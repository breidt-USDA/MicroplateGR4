function stats = LinReg(XYdata)
% Function: stats = LogReductLinear(LogCFUdata)
% Version 08132010
% Purpose: Linear regression of XY data. The algorithm was written by: 
% Dr. Fred Breidt, USDA/ARS
% microbiologist at NC State University
% Dr. Jason A. Osborne, 
% NC State, University Department of Statistics, Raleigh, NC. 
% 
% Contact information: 
% Fred Breidt
% USDA/ARS Microbiologist
% Department of Food Bioprocessing and Nutrition Sciences
% 322 Schaub Hall, Box 7624
% NC State University, Raleigh, NC 27695-7624
% TEL: (919) 513-0186
% 
% Function Input Parameters
% XYdata = an N x 2 vector with the first column of data being X values 
%   and the second column Y values (typically time, Ln(OD)) 
% 
% Return values: stats data structure
% stats.Time                    1 x n array of time values
% stats.LogCFUpermL             1 x n array of CFU/ml (survivors)
% stats.PredCFUpermL            1 x n array of predicted CFU/ml values
% stats.Slope                   slope of the regression line
% stats.StdErrSlope             standard error for slope estimate
% stats.Intercept               intercept (initial CFU/ml)
% stats.StdErrIntercept         std. error for intercept estimate
% stats.Residuals               error values for each data point
% stats.Rsquared                regression coefficient (R squared value)
%--------------------------------------------------------------------------

% setup the design and Y matrices. The input data structure, 
% XYdata, has two columns wich are seperated into the design matrix (first 
% column is ones, second column is X data) and a column vector of Y data

datapoints = size(XYdata);              %get size of input matrix
Npts = datapoints(1);                   %number of data points
X = ones(Npts,2);                       %initialize n x 2 design matrix 
X(:,2) = XYdata(:,1);                   %load time data from input matrix
Y = XYdata(:,2);                        %load Ydata

% get predicted regression parameters (beta hat, 2 x 1 matrix). Bhat(1) is
% the Y intercept and Bhat(2) is the slope.
Bhat = ((X'*X)^-1)*X'*Y;            

%determine the mean square error. 
MSE =  ((Y-(X*Bhat))'*(Y - (X*Bhat)))/(Npts - 2);

%determine the variance-covariance matrix (VCV, a 2 x 2 matrix) with the 
% variance of the parameters on the diagonal: (1,1) is intercept variance,
% (2,2) is slope varience.
VCV = MSE*((X'*X)^-1);

%determine predicted values.
pred = X*Bhat;                          %design matrix times param vector

%determine the residuals
resid = Y - pred;                       %Y values - predicted

%coefficient of determination (Rsquared value) 
Ybar = sum(Y)/Npts;                     %mean of Y values
SSM = sum((pred-Ybar).*(pred-Ybar));    %predicted minus mean squared
SSE = sum((Y-pred).*(Y-pred));          %Y vals minus pred squared (error)
SST = SSM+SSE;                          %sum squares total
Rsq = SSM/SST;                          %r squared value

%set up return structure
stats.X = X(:,2);                       %array of X values
stats.Y = Y;                            %array of Y values
stats.PredY = pred;                     %array of Predicted Y for X
stats.Slope = Bhat(2);                  %slope of the regression line
stats.StdErrSlope = sqrt(VCV(2,2));     %standard error for slope estimate
stats.Intercept = Bhat(1);              %intercept
stats.StdErrIntercept = sqrt(VCV(1,1)); %std. error for intercept estimate
stats.Residuals = resid;                %residuals
stats.Rsquared = Rsq;                   %r squared value
end

