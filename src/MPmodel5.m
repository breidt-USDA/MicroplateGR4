classdef MPmodel5 < fitGRModel
   % F. Breidt USDA/ARS MPmodel5 Version 5.0
   % Processing optical density readings from a 96 well microtiter plate to
   % generate microbial growth kinetics data 
   %
   % MPmodel5 contains the core functions for converting raw microtiter
   % plate data files to an organized 3D matrix structure (8 x 12 x N), 
   % for 8 rows, 12 columns and N timepoints. Raw data consists of a Matlab
   % table (N x 96) obtained with time - OD data for each well in
   % row-dominant format (A1, A2, A3, ... H12). Growth kinetics data are 
   % generated from one or more user-selected groups of wells (for both 
   % blanks and data wells) using SW ("sliding window") processing method.
   %
   % References:
   % 1. Breidt, F., T. L. Romick, and H. P. Fleming. 1994. A rapid method 
   % for the determination of bacterial growth kinetics. J. Rapid Methods 
   % Automation Microbiol. 3:59-68.
   % https://doi.org/10.1111/j.1745-4581.1994.tb00308.x
   % 2. Atolia, E., Cesar, S., Arjes, H. A., Rajendram, M., Shi, H., Knapp, 
   % B. D., Khare, S., Aranda-Diaz, A., Lenski, R. E., & Huang, K. C. 
   % (2020). Environmental and physiological factors affecting high 
   % throughput measurements of bacterial growth. mBio, 11, 
   % e01378-20. https://doi.org/10.1128/mBio.01378-20 .  
   %----------------------------------------------------------------------
   properties   
      MP_XYY      %Microplate data matrix 
      MP_XYN      %Microplate in 8x12xN data matrix format 
      TPs         %number of time points (including time zero)
      MachineRes  %resolution of spectrophotometer
      FirstN      %initial points (1..N) for determining min OD val
      blankData   %info about the blank value used for subtracting from OD
      fitmodel    %model fitting class fitGRModel object
      default_model %Default curve fitting: Gompertz (1) or Baranyi (2) 
      dataname    %string with name of data object
      MPmodel5_version  %1.1 fixed PERCENT to = cormin*100/uncormin
   end % properties

   methods %public methods------------------------------------------------
      
      %constructor with microplate data and text string (for data name)
      function obj = MPmodel5(MP_XYYtable, dataname)
     
         %take Nx97 Table and make MP_XYY matrix, MP_XYN 3D matrix 
         obj@fitGRModel(); %curve fitting
         obj.default_model = 1; %Gompertz (1),Baranyi (2)
         obj.MP_XYY = MP_XYYtable{:,:};      
         obj.TPs = size(obj.MP_XYY,1);       
         obj.MP_XYN = zeros(8,12,obj.TPs);   
         obj.setXYN(0);                      
         obj.MachineRes = 0.001; %based on spectrophotometer resolution
         obj.FirstN = 3;                            
         obj.blankData = struct([]); 
         obj.dataname = dataname;            
         obj.MPmodel5_version = 1.1;             
      end %function

      %sets the growth model to use for graphing and for optimization 
      function setModel(obj,value)
         
         obj.default_model = value;
      end

      %Assume plateinfo is filled
      function blankdata = getBlankData(obj, blankCellStr)
         
         %get blank data from a cell array with well identifiers
         if ~strcmp(blankCellStr,"")
            blankmat = obj.convertCellStr(blankCellStr);
            blankXYY = obj.getXYYMat(blankmat); 
            XY = obj.XYY2XY(blankXYY);
            blankdata.XY = XY;
            blankdata.XYY = blankXYY;
            blankdata.labelstr = blankCellStr;
            blankdata.Nvals = size(XY,1);
            blankdata.mean = mean(XY(:,2));
            blankdata.median = median(XY(:,2));
            blankdata.mode = mode(XY(:,2));
            blankdata.std = std(XY(:,2));
         else
            blankdata.XY = [0,0];
            blankdata.XYY= [0,0];
            blankdata.labelstr = "";
            blankdata.Nvals = 1;
            blankdata.mean = 0;
            blankdata.median = 0;
            blankdata.mode = 0;
            blankdata.dtd = 0;
         end
         obj.blankData = blankdata;
      end %function

      %get a table of results for each curve in rowcol data matrix
      function ResTable = getResultsTable(obj, XYY, LnXYY, ...
            curveCellStr, blankmean, processbyN)
         
         %build matrix with results for each row and column
         % GRdata.ResTable has: 'GR','DBL','LAG','MIN','MAX' columns;
         Xdata = LnXYY(:,1);
         Nrc = size(LnXYY,2); 
         vNames = cellstr(["GR","DBL","LAG","MIN","MAX"]); 
         vTypes = cellstr(["double","double","double","double","double"]);
         sz = [Nrc-1,5]; % subtract the time column
         TempTable = table('Size',sz,'VariableNames',vNames,...
            'VariableTypes', vTypes);
         CorMinTab = table('Size',[Nrc-1,1],'VariableNames',"COR_MIN", ...
             'VariableTypes',cellstr("double"));
         for i=2:Nrc %get data for each curve                         
            temp = obj.getGRdata([Xdata, LnXYY(:,i)], processbyN);
            TempTable(i-1,:) = temp.ResTable;
            % reset min and max to corrected real values (minus blank)
            TempTable.MIN(i-1) = min(XYY(:,i));
            CorMinTab.COR_MIN(i-1) = min(XYY(:,i)) - blankmean; 
            TempTable.MAX(i-1) = max(XYY(:,i))-blankmean;
            %TempTable.RMSE(i-1) = temp.RMSE;
         end %for
         %Add row and column values to the table       
         ResTable = [curveCellStr', TempTable];      %add to data
         ResTable.Properties.VariableNames(1) = {'RowCol'};
         ResTable = renamevars(ResTable,'MIN','UNCOR_MIN');
         ResTable = renamevars(ResTable,'MAX','COR_MAX');
         %corrected min OD as a percent of the blank value
         if blankmean ~= 0
            PERCENT = CorMinTab.COR_MIN*100./ResTable.UNCOR_MIN;
         else
            PERCENT = CorMinTab.COR_MIN*0; %vector of zeros
         end
         ResTable.PERCENT = PERCENT;
         ResTable = [ResTable,CorMinTab];
         GRmean = mean(ResTable.GR); 
         GRstd = std(ResTable.GR);
         DBLmean = mean(ResTable.DBL);
         DBLstd = std(ResTable.DBL);
         LAGmean = mean(ResTable.LAG);
         LAGstd = std(ResTable.LAG);
         MINmean = mean(ResTable.UNCOR_MIN);
         MINstd = std(ResTable.UNCOR_MIN);
         MAXmean = mean(ResTable.COR_MAX);
         MAXstd = std(ResTable.COR_MAX);
         PERCENTmean = mean(ResTable.PERCENT);
         PERCENTstd = std(ResTable.PERCENT);
         CMINmean = mean(CorMinTab.COR_MIN);
         CMINstd = std(CorMinTab.COR_MIN);
         %add two rows with means and one with model from all data
         temp = obj.getGRdata(LnXYY,processbyN); %one col many vals
         allres = temp.ResTable;
         newrows = array2table([GRmean,DBLmean,LAGmean,MINmean,MAXmean, ...
            PERCENTmean,CMINmean; ... 
            GRstd,DBLstd,LAGstd,MINstd,MAXstd,PERCENTstd,CMINstd; 
            allres.GR,allres.DBL,allres.LAG,NaN, ...
            exp(allres.MAX),NaN,exp(allres.MIN)]);
         newrows.Properties.VariableNames = {'GR','DBL','LAG','UNCOR_MIN', ...
            'COR_MAX','PERCENT','COR_MIN'};
         nametab = cell2table({'mean';'std';'all'},'VariableNames', ...
            {'RowCol'});
         newrows = [nametab, newrows];
         ResTable = [ResTable;newrows];
         ResTable = movevars(ResTable,"COR_MIN",'After',"UNCOR_MIN");
      end %function       

   end % end of public functions ------------------------------------------

   methods (Access = protected)

      %convert a cellstring vector to data matrix 2xN
      function datamat = convertCellStr(obj,cellStrVec)

         %for each cellStr value, convert alphanumeric to number matrix
         if ~isempty(cellStrVec)
            vecLen = length(cellStrVec);
            datamat = zeros(2,vecLen);
            for i=1:vecLen
               tempstr = cellStrVec{i};
               datamat(1,i) = obj.getNumberVal(tempstr(1)); % row letter
               datamat(2,i) = str2double(tempstr(2:end));
            end % for
         else
            datamat = [];
         end

      end % function
    
      %Return structure with growth rate params and stats
      function GRdata = getGRdata(obj, LnXYY, processbyN)
      
         %uses CORRECTED natural log OD data values
         [Ndp, Ncols] = size(LnXYY);         
         LnXY = obj.XYY2XY(LnXYY);          
         Nys = Ncols - 1;                    %number of curves in data set
         istart = 1;                        
         iend = processbyN * Nys;            %initial 'window' for slope
         LR = obj.getLRegdata(LnXY,istart,iend); 
         startindex = istart;                %save initial start index
         istart = istart + Nys;              
         iend = iend + Nys;                  
         %for each window through the growth curve calculate the slope
         while iend <= Ndp                     
            LRtemp = obj.getLRegdata(LnXY,istart,iend); 
            if LRtemp.Slope > LR.Slope       %check for largest slope
               startindex = istart;          
               LR = LRtemp;                  %save linear regression data            
            end %if
            istart = istart + Nys;           
            iend = iend + Nys;               
         end %while loop
         %find the minimum slope, and minimum and maximum OD
         minOD = mean(LnXY(1:Nys,2));        %mean of initial values  
         maxOD = obj.getMaxOD(LnXYY,Nys);    %mean of Max from each curve         
         DblTime = 0;                        
         LagTime = 0;                        
         if LR.Slope > 0                     %check for positive slope
            DblTime = log(2)/LR.Slope;       
            %if minOD > log(obj.MachineRes)
               LagTime = (minOD - LR.Intercept)/LR.Slope; %calc lag time
            %end
         end %if 
         %calculate XY values for the model from sliding window paraeters
         pvec = [minOD,LR.Slope,LagTime,maxOD];
         Xvals = LnXY(:,1);
         if obj.default_model == 1
            modelXY = obj.calcGompertz(Xvals,pvec);
            model_name = "Gompertz";
         else
            modelXY = obj.calcBaranyi(Xvals,pvec);
            model_name = "Baranyi";
         end %if
         %fill the GR data structure with regression data
         GRdata.Model = model_name;
         GRdata.ObsLnXYY = LnXYY;
         GRdata.ObsLnXY = LnXY;
         GRdata.PredLnXY = modelXY;
         DataVec = [LR.Slope,DblTime,LagTime,minOD,maxOD];    
         GRdata.ResTable = array2table(DataVec); 
         GRdata.ResTable.Properties.VariableNames = {'GR','DBL','LAG',...
            'MIN','MAX'};                    
         GRdata.Npts = Ndp;                  
         GRdata.IndexStart = startindex;  %start and end for max GR window 
         GRdata.IndexEnd = startindex + processbyN -1; 
         GRdata.MinOD = minOD;
         GRdata.MaxOD = maxOD;
         GRdata.MachineResolution = obj.MachineRes;
         GRdata.LogMachineResolution = log(obj.MachineRes);
         GRdata.Slope = LR.Slope;            
         GRdata.YIntercept = LR.Intercept;   
         GRdata.Rsq = LR.Rsquared;           
         GRdata.Residuals = LR.Residuals;
         fit = obj.getFit(LnXY,GRdata.PredLnXY,Ndp);
         GRdata.SSE = fit.SSE;
         GRdata.RMSE = fit.RMSE;
      end %function

      %get Sum Squared Error and Root Mean Square Error
      function fit = getFit(~, observed, predicted, Ndp)
         
         %statistics for fitted data
         fit = struct;
         error = 0;
         if Ndp > 0
            for i=1:Ndp
               error = error + (observed(i,2) - predicted(i,2))^2;
            end
            fit.SSE = error;
            fit.RMSE = sqrt(error/Ndp);
         end %if
      end %function
        
      %used for output of writeRowColTable
      function updatedTable = Row2Letters(obj, restable)
         
         %convert first col of restable from numbers to letters  
         numbercol = restable{:,1};
         DATA = string(numbercol);
         Nrc = height(restable);
         for j=1:Nrc
            DATA(j) = obj.getRowLetter(DATA(j));
            DATA(j) = strcat(DATA(j),string(restable{j,2}));
         end %for
         temptab = array2table(DATA);
         updatedTable = [temptab, restable(:,3:end)];
      end % end
  
   
      %MODIFIED: convert all submatrix OD values to corrected natural log
      function LnXYY = getLnData(obj,XYYmat,blankval)
      
         %subtract the blank (replace negative or zero) then take natl log 
         submat = XYYmat(:,2:end) - blankval;
         submat = max(submat,obj.MachineRes); %REPLACE all below res
         XYYmat(:,2:end) = log(submat);
         LnXYY = XYYmat;
      end %function

      %do linear regression for sub set of data (start to end X vals)
      function LRegdata = getLRegdata(~,LnXY,starttime,endtime)
      
         % make XY matrix for linear regression (LinReg)
         dataXvals = LnXY(starttime:endtime,1);  
         dataYvals = LnXY(starttime:endtime,2);  
         datamat = [dataXvals dataYvals];    
         LRegdata = LinReg(datamat); 
      end %function

      %get max of all OD values
      function maxOD = getMaxOD(~,XYY,Nys)
         
         maxvec = zeros(1,Nys);
         for i=1:Nys
            maxvec(i) = max(XYY(:,i+1));
         end
         maxOD = mean(maxvec);
      end %function

      % NEW get XYY from a data mat 2xN matrix      
      function XYYmat = getXYYMat(obj,datamat)

         %create and fill XYY
         Nwells = size(datamat,2);
         XYYmat = zeros(obj.TPs,Nwells+1);
         XYYmat(:,1) = obj.MP_XYY(:,1); 
         for i=1:Nwells
            tempXY = obj.getXY(datamat(1,i),datamat(2,i));
            XYYmat(:,i+1) = tempXY(:,2);
         end %for
      end %function

      % USED get curve from a perticular row and col in the 8 x 12 microplate
      function XY = getXY(obj,row,col)
      
         %extract row-col as Nx2 from 3D XYN matrix
         yvals = squeeze(obj.MP_XYN(row,col,1:end)); 
         xvals = obj.MP_XYY(:,1);           
         XY = [xvals yvals];                 
      end %function

      %get one merged Nx2 from matrix of X with multiple Y values 
      function XYall = XYY2XY(~,XYY)
          
         %loop through Y columns to construct single Nx2
         [Nrows, Ncols] = size(XYY);         
         NYcols = Ncols - 1;
         XYlen = Nrows*NYcols;              
         resindex = 1;
         XYall = zeros(XYlen,2);
         %manually rebuild data into a single Nx2
         for i=1:NYcols                      
            for j=1:Nrows                    
               XYall(resindex,1) = XYY(j,1);
               XYall(resindex,2) = XYY(j,i+1); 
               resindex = resindex + 1;      
            end %for
         end %for
         XYall = sortrows(XYall,1);
      end %function

      %build an XYN (row, col, OD value) 3-D matrix from XYY 2-D matrix
      function setXYN(obj, blankval)
         
         %loop through XYY and for each row, col add OD vector
         colindex = 2;                       
         for i=1:8                           
            for j=1:12                       
               obj.MP_XYN(i,j,:) = obj.MP_XYY(:,colindex)-blankval; 
               colindex = colindex + 1;      
            end %for
         end %for 
      end %function
 
      %get text string for graph labels with rows and columns used
      function Xlab = getLabelVec(obj,rowcolmat)

         %make labels from a row and column matrix
         Xlabels = num2str(rowcolmat);
         Nvals = size(Xlabels,1);
         for i=1:Nvals
            Xlabels(i,1) = obj.getRowLetter(Xlabels(i,1)); 
         end
         Xlab = erase(string(Xlabels), " ");  %remove spaces from label
      end

      %convert a numeric well number to a row-col cell string
      function cellStr = getRowColStr(obj,wellNumberVec)
         
         %well number = 1 - 96, rowcol strings = {'A1','B3', etc...}
         nvals = length(wellNumberVec);
         cellStr = cell(1,nvals);
         for i=1:nvals
            val = wellNumberVec(i);
            cellstrval = obj.well2RowCol(val);
            cellStr(i) = cellstr(cellstrval);
         end
      end
      
      %convert a well number to a row-col string
      function RowColStr = well2RowCol(~,well_number)

          % well numbers are linear from 1:96 (across)
          col = mod(well_number, 12);
          if col == 0
            col = 12;
          end
          colstr = string(col);
          checkval = well_number/12;
          if checkval <=1
              RowColStr = "A" + colstr;
              return;
          elseif checkval <= 2
              RowColStr = "B" + colstr;
              return;
          elseif checkval <= 3
              RowColStr = "C" + colstr;
              return;
          elseif checkval <= 4
                  RowColStr = "D" + colstr;
                  return
          elseif checkval <= 5
                  RowColStr = "E" + colstr;
                  return;
          elseif  checkval <= 6
                  RowColStr = "F" + colstr;
                  return;
          elseif  checkval <= 7
                  RowColStr = "G" + colstr;
                  return;
          else 
             RowColStr = "H" + colstr;
          end %if
      end %function
  
      %convert number char to a capital letter char
      function Xchar = getRowLetter(~,charval)  

         switch charval
            case "1"
               Xchar = 'A';
            case "2"
               Xchar = 'B';
            case "3"
               Xchar = 'C';
            case "4"
               Xchar = 'D';
            case "5"
               Xchar = 'E';
            case "6"
               Xchar = 'F';
            case "7"
               Xchar = 'G';
            case "8"
               Xchar = 'H';
         end %switch
      end %function

      %conert letter to corresponding number
      function value = getNumberVal(~,strletter)

         switch (strletter)
            case 'A' 
               value = 1;
            case 'B'
               value = 2;
            case 'C'
               value = 3;
            case 'D'
               value = 4;
            case 'E'
               value = 5;
            case 'F'
               value = 6;
            case 'G'
               value = 7;
            case 'H'
               value = 8;
            otherwise
               value = 9; %not in A->H
         end % switch
      end % function

   end % private methods -----------------------------------------------

end %class def