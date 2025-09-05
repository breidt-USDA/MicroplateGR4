classdef ProcWells2 < MPmodel5
   % F. Breidt USDA/ARS ProcWells Version 1.0
   % ProcWells has functions for processing a selected set of wells in a
   % microtiterplate, processing and displaying OD data and growth kinetics
   % MPmodel5 is the parent class with plate processing functions
   % fitGRModel is base clase (optimization) < handle class
   %----------------------------------------------------------------------
   properties 
      wellList
      ProcPlate_Version
   end % properties

   methods

      %constructor
      function obj = ProcWells2(MP_XYYtable, dataname)

         %call the base class constructor, set version number
         obj@MPmodel5(MP_XYYtable, dataname);
         obj.ProcPlate_Version = 1.0;
         obj.wellList = zeros(96,1);
      end % function

      % Get tile number based on the axes selected by mouse click
      function setWellList(obj, t, ax)
        
         % t is the tile graph, ax = a selected graph in t
         tileIndex = find(t.Children == ax);
         if sum(ax.Color - [1,0.5, 0.5]) == 0
            set(ax, 'Color', [1,1,1]);
            obj.wellList(tileIndex) = 0;
         else
            set(ax, 'Color', [1, 0.5, 0.5]);
            obj.wellList(tileIndex) = 1;
         end %if 
      end %function

      % remove selections and clear the selection list
      function clearWells(obj)

         obj.wellList(obj.wellList ~= 0) = 0; %zero wellList
         obj.showPlate();
      end % function

      %get a list string from selections on graph
      function listString = getListString(obj)
         
         %get cell array with row-col data e.g. {'A1'} {'C2'}
         nwells = sum(obj.wellList);
         listString = "";
         index = 1;
         if nwells > 0
            tempvec = zeros(nwells,1);
            for i = 1:96
               if obj.wellList(i) > 0
                  tempvec(index) = i;
                  index = index + 1;
               end % if
            end % if
            listString = obj.getRowColStr(tempvec);
         end % if
      end % function

      %get LnXY data for regression line used for growth rate
      function slopeline = getSlopeLine(~,minLn,maxLn,slope,Yintercept)
         
         %returns a Nx2 line for plotting based on input parameters
         lineY = linspace(minLn,maxLn); 
         if slope ~= 0
            lineX = (lineY - Yintercept)./slope;
            slopeline = [lineX',lineY'];
         else
            slopeline = [0,0];
         end %if
      end % function

      %check PDstruct for valid parameters flag
      function tf = validPlateInfo(obj)
         
         %return true-false value for valid plate data
         tf = obj.PDstruct.valid;
      end % function

      % return CrvData  and GR data from one growth curve
      function CrvData = getCurveData(obj, curveCellStr, blankCellStr, ...
            processbyN)

         %get gr data and results table
         CrvData.blankdata = obj.getBlankData(blankCellStr);
         CrvData.blankCellStr = strjoin(blankCellStr,',');
         CrvData.dataCellStr = strjoin(curveCellStr,',');
         curvemat = obj.convertCellStr(curveCellStr);
         XYYmat = obj.getXYYMat(curvemat);
         CrvData.XYY = XYYmat;
         CrvData.XY = obj.XYY2XY(XYYmat);
         LnXYY = obj.getLnData(XYYmat,CrvData.blankdata.mean);
         CrvData.LnXYY = LnXYY; %data are CORRECTED for blank
         CrvData.LnXY = obj.XYY2XY(LnXYY);
         CrvData.grdata = obj.getGRdata(LnXYY,processbyN);
         CrvData.ResultsTable = obj.getResultsTable(XYYmat, LnXYY, ...
            curveCellStr, CrvData.blankdata.mean,processbyN);
      end %function

      function [midtime, midval] = getMidPoint(~,curveData)
         startindx = curveData.grdata.IndexStart;
         endindx = curveData.grdata.IndexEnd;
         starttime = curveData.LnXY(startindx,1);
         endtime = curveData.LnXY(endindx,1);
         startval = curveData.LnXY(startindx,2);
         endval = curveData.LnXY(endindx,2);
         midtime = (starttime + endtime)/2;
         midval = (startval + endval)/2;      
      end

      %show microplate, for diagnostic purposes 
      function tilegraph = showPlate(obj)

         %set up at figure
         obj.wellList(obj.wellList ~= 0) = 0; %zero wellList
         t = tiledlayout(8,12,'TileSpacing','none'); %8x12 grid, no spacing   
         xlabel(t,"Time (h)",'FontSize',18,'FontWeight','bold');              
         ylabel(t,"Optical Density (600 nm)",'FontSize',18,'FontWeight', ...
            'bold'); 
         Ydata = obj.MP_XYY(:,2:end);        % just get Y data values
         ymax_val = max(Ydata,[],"all");     % max Y val for Y axis scale
         for i=2:97                          
            ax = nexttile(t);                         
            plot(ax,obj.MP_XYY(:,1),obj.MP_XYY(:,i),'.k'); % plot in tile
            ylim(ax, [0,ymax_val]);              
            set(ax,'xtick',[],'ytick',[]);
            %ax.Toolbar.Visible = 'off'; %very slow execution with 2025A
            ax.Toolbar = []; % MATLAB Fix
            ax.ButtonDownFcn = @(src, event) obj.setWellList(t,ax);
         end %for
         rowvals = ['A','B','C','D','E','F','G','H'];
         colvals = ["1","2","3","4","5","6","7","8","9","10","11","12"];
         t.Children = flipud(t.Children);
         for j=1:8
            rnum = tilenum(t,j,1);
            ylabel(t.Children(rnum),rowvals(j),'Rotation',0, ...
               'FontSize',14,'FontWeight','bold');
         end %for
         for k=1:12
            cnum = tilenum(t,1,k);
            title(t.Children(cnum),colvals(k),'FontSize',14, ...
               'FontWeight','bold');
         end %for
         title(t,obj.dataname,'FontSize',18,'FontWeight','bold'); %title
         tilegraph = t;
      end %function

      function plotBlankData(~,blankdata)

         % Histogram for blank data values
         fh = figure;
         ax = axes('Parent',fh);
         histogram(blankdata.XY(:,2)); %plot histogram
         blankval = sprintf('%.3g',blankdata.mean);
         stdval = sprintf('%.2g',blankdata.std);
         hist_titlestr = append('Histogram, ', ...
         strjoin(blankdata.labelstr,","), ...
         newline, ' mean (stdev): ', ...
         blankval,' (',stdval, ')');
         title(hist_titlestr, 'FontSize',18);
         ax.XAxis.FontSize = 14;
         ax.YAxis.FontSize = 14;
         ax.XAxis.FontWeight = 'bold';
         ax.YAxis.FontWeight = 'bold';
         ax.XAxis.TickLabelFormat = '%.3f';
         ax.YAxis.TickLabelFormat = '%.0f';
         %set label font only
         xlabel('OD value','FontSize',18, 'FontWeight', 'bold');
         ylabel('Number', 'FontSize',18, 'FontWeight', 'bold');
      end % function
   
      %get Ln XY graph with observed and predicted data
      function crvdata = showModel(obj,curveCellStr,blankCellStr, ...
            processbyN)
         
         %main graph of SW model data with slope, min, max lines
         crvdata = obj.getCurveData(curveCellStr,blankCellStr, ...
            processbyN);
         minY = crvdata.grdata.ResTable.MIN;
         maxY = crvdata.grdata.ResTable.MAX;
         slope = crvdata.grdata.Slope;
         Yintercept = crvdata.grdata.YIntercept;
         sline = obj.getSlopeLine(minY,maxY,slope,Yintercept);
         XYdata = crvdata.grdata.ObsLnXY;
         XYpred = crvdata.grdata.PredLnXY;
         Ndp = size(XYpred,1);
         if Ndp > 0
            SSE = sum((XYdata(:,2)-XYpred(:,2)).^2);
            val = sqrt(SSE/Ndp);
            RMSE = sprintf('%.3g',val);
         else
            RMSE = "nd";
         end 
         fh = figure;
         ax = axes('Parent',fh);
         hold on
         plot(sline(:,1),sline(:,2),'-k',"LineWidth",2)
         plot(XYdata(:,1),XYdata(:,2),'.k', ...
            XYpred(:,1),XYpred(:,2),'--r', 'LineWidth',2);
         yline([minY, maxY],'--b','LineWidth',2);
         titlestr = strjoin(curveCellStr,",");
         mname = crvdata.grdata.Model;
         titlestr = titlestr + ": " + mname + " (SW), RMSE = " + RMSE;
         title(titlestr,'FontSize',18);
         ax.XAxis.FontSize = 14;
         ax.YAxis.FontSize = 14;
         ax.XAxis.FontWeight = 'bold';
         ax.YAxis.FontWeight = 'bold';
         xlabel('Time (h)','FontSize',18,'FontWeight','bold');
         ylabel('Ln Optical Density','FontSize',18,'FontWeight','bold');
         xlim([0,XYdata(end,1)]);
         ax.XAxis.TickLabelFormat = '%.0f';
         ax.YAxis.TickLabelFormat = '%.1f';
         ax.YAxis.MinorTick = 'on';
         box on
         hold off
      end % function

   
      %plot observed and predicted data (predCD = curveData struct)
      function plotfit(obj,XYdata,XYpred,pred_params,predCD, ...
            graphtitle, RMSE)

         %params: MinOD,GR,lag,MaxOD
         minY = pred_params(1);
         slope = pred_params(2);
         lag = pred_params(3);
         Yintercept = minY - (slope*lag);
         maxY = pred_params(4);
         sline = obj.getSlopeLine(minY,maxY,slope,Yintercept);
         RMSEstr = sprintf('%.3g',RMSE);
         fh = figure;
         ax = axes('Parent',fh);
         hold on
         plot(sline(:,1),sline(:,2),'-k', "LineWidth",2)
         plot(XYdata(:,1),XYdata(:,2),'.k',XYpred(:,1),XYpred(:,2), ...
            '--r','LineWidth',2);
         yline([minY, maxY],'--b','LineWidth',2);
         mname = predCD.grdata.Model;
         titlestr = graphtitle + ": " + mname + " (fit),";
         titlestr = titlestr + " RMSE = " + RMSEstr;
         title(titlestr, 'FontSize',18);
         ax.XAxis.FontSize = 14;
         ax.YAxis.FontSize = 14;
         ax.XAxis.FontWeight = 'bold';
         ax.YAxis.FontWeight = 'bold';
         ax.YAxis.MinorTick = 'on';
         xlabel('Time (h)', 'FontSize',18,'FontWeight','bold');
         ylabel('Ln Optical Density', 'FontSize',18,'FontWeight', 'bold');
         xlim([0,XYdata(end,1)]);
         box on
         hold off
      end

      %show graph of selected curves using corrected XY OD data (not Ln)
      function showCurves(obj,curveCellStr)

         %observed data, OD curves of selected wells (not corrected)
         fh = figure;
         ax = axes('Parent',fh);
         curvemat = obj.convertCellStr(curveCellStr);
         XYYmat = obj.getXYYMat(curvemat);
         XYdata = obj.XYY2XY(XYYmat);
         plot(XYdata(:,1),XYdata(:,2),'ok');
         titlestr = strjoin(curveCellStr,",");
         titlestr = titlestr + ": OD Data";
         title(titlestr, 'FontSize',20);
         ax.XAxis.FontSize = 14;
         ax.YAxis.FontSize = 14;
         ax.XAxis.FontWeight = 'bold';
         ax.YAxis.FontWeight = 'bold';
         ax.YAxis.MinorTick = 'on';
         xlabel('Time (h)', 'FontSize',18,'FontWeight','bold');
         ylabel('Optical Density', 'FontSize',18,'FontWeight', 'bold');
         box on;
      end % function
  
   end % public methods

end %end class
