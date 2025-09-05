classdef fitGRModel < handle
%  FB fitGRModel Version 1.0, Matlab base class (handle class)
%  Constructor takes no arguments, separate functions fit observed data 
%   matrix with time, natural log OD to either Gompertz or Baranyi model 
%  Each function requires OD data and initial parameter vector
%  Separate functions for calculation of model curves given parameters
% 
% MODEL REFERENCES: 
%  Baranyi: Grijspeerdt and Vanrolleghem, 1999, Food Microbiol 16,593-605
%  Gompertz: Zwietering et al, 1990, Appl Env Microbiol 56(6),1875-1881
%  NOTE: To get statistics on estimated parameters a bootstrapping method
%  may be used - currently not implemented. 
%-------------------------------------------------------------------------

   properties
      obsData     % Nx2, time and OD
      predData    % Nx2, time and OD
   end %properties

   methods 

      % Constructor, call setModel with empty matrix
      function obj = fitGRModel()
         
         obj.obsData = [];
         obj.predData = [];
      end %function

      % Gompertz model takes init params, returns opt params
      function res = optimizeGompertz(obj,obsData,initParams)

         %set initial obs and pred data
         obj.obsData = obsData;
         fh = @obj.calcGompertzSSE; %returns SSE uses Obs XY, Pred data
         [params,Fval,exitFlag,outPut] = fminsearch(fh, initParams,  ...
            optimset('PlotFcn', @optimplotfval, 'TolX',1e-6,'MaxFunEvals', ...
            1e6,'MaxIter', 1e6,'TolFun', 1e-6));
         res = obj.setResults(initParams,params,Fval,exitFlag,outPut);
      end % function

      % Baranyi model takes init params, returns opt params
      function res = optimizeBaranyi(obj, obsData,initParams)
         
         %fminsearch optimization using Baranyi model
         %paramvec: initVal, GR, lagTime, maxOD
         obj.obsData = obsData;
         fh = @obj.calcBaranyiSSE;
         [params,Fval,exitFlag,outPut] = fminsearch(fh, initParams,  ...
            optimset('PlotFcn', @optimplotfval,'TolX',1e-6,'MaxFunEvals', ...
            1e6,'MaxIter', 1e6, 'TolFun', 1e-6));
         res = obj.setResults(initParams,params,Fval,exitFlag,outPut);
      end % function 

      %from Grijspeerdt and Vanrolleghem 1999 Food Microbiol. 16,593-605
      function XYdata = calcBaranyi(~, timeVec, params)
         
         %calculate predicted XY data for the Baranyi model
         %paramvec: initVal, GR, lagTime, maxOD
         Yo = params(1);
         Umax = params(2);
         lambda = params(3);   %note: lag (lambda) = ho/Umax
         Ymax = params(4);
         Ndp = length(timeVec);              %number of timepoints
         predicted = zeros(Ndp,2);             
         predicted(:,1) = timeVec;                                 
         predicted(1,2) = Yo;               %assign initial Y value
         for i=2:Ndp 
            At = timeVec(i) + (1/Umax)*log(exp(-1*Umax*timeVec(i))+ ...
               exp(-1*Umax*lambda)-exp(-1*Umax*(timeVec(i)+lambda)));
            predicted(i,2) = Yo + Umax*At  - log(1 + ...
               (exp(Umax*At)-1)/(exp(Ymax - Yo))); 
         end %for 
         XYdata = predicted;
      end %function

      %From Zwittering et al 1992 paper, JFP
      function XYdata = calcGompertz(~,timeVec,params)
         
         %calculate Gompertz model to generate XY data from parameters
         %paramvec: initVal, GR, lagTime, maxOD
         Yo = params(1);
         Umax = params(2);
         lambda = params(3);
         Ymax = params(4);
         Ndp = length(timeVec);              %number of timepoints
         predicted = zeros(Ndp,2);             
         predicted(:,1) = timeVec;                                
         A = Ymax + abs(Yo); 
         e = exp(1);                                         
         mu = Umax;                    
         predicted(1,2) = Yo;               %assign initial Y value
         %calculate perdicted Y value for each timepoint in the curve
         for i=2:Ndp                         
            predicted(i,2) = Yo + ...       
               A*exp(-1*exp((mu*e/A)*(lambda - timeVec(i))+1));
         end %for 
         XYdata = predicted;
      end %function
  
      %plot observed and predicted data
      function plotObsPred(obj)
         
         %simple observed and predicted data plot
         plot(obj.obsData(:,1),obj.obsData(:,2),'ok', ...
            obj.predData(:,1),obj.predData(:,2),'-r');
      end %function

   end % end public methods

   methods (Access = private) %-------------------------------------------

      %from Grijspeerdt and Vanrolleghem 1999 Food Microbiol. 16,593-605
      function SSE = calcBaranyiSSE(obj, params)
         
         %return SSE for Baranyi fit of observed data (obj.ObsData)
         timevec = obj.obsData(:,1);
         obj.predData = obj.calcBaranyi(timevec,params);
         SSE = obj.getSSE;
      end %function

      %From Zwittering et al 1992 paper, JFP
      function SSE = calcGompertzSSE(obj,params)
         
         %return SSE for Gompertz model with observed (obj.ObsData) 
         timevec = obj.obsData(:,1);
         obj.predData = obj.calcGompertz(timevec,params);
         SSE = obj.getSSE;
      end %function

      % calc sum of squared error term
      function SSE = getSSE(obj)

         temp = (obj.obsData(:,2) - obj.predData(:,2)).^2;
         SSE = sum(temp);
      end % function
   
      % calc root mean squared error
      function RMSE = getRMSE(obj)
         
         ndp = size(obj.obsData,1);
         SSE = obj.getSSE();
         if ndp > 0
            RMSE = sqrt(SSE/ndp);
         else
            RMSE = 0;
         end % if
      end % function

      % set up return values in data structure
      function res = setResults(obj,initParams,params,Fval,exitFlag,outPut)

         %a structure with results from fminsearch
         res.obsData = obj.obsData;
         res.predData = obj.predData;
         res.initParams = initParams;
         res.predParams = params;
         res.SSE = Fval;
         res.RSME = obj.getRMSE();
         res.exitFlag = exitFlag;
         res.funcCount = outPut.funcCount;
         res.message = outPut.message; 
      end % function
      
   end % end private methods

end %classdef