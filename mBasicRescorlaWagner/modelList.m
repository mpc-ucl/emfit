function model = modelList; 
% 
% A file like this should be contained in each model class folder and list the
% models to be run, together with some descriptive features. 
% 
% GENERAL INFO: 
% 
% list the models to run here. The models must be defined as likelihood functions in
% the models folder. They must have the form: 
% 
%    [l,dl,dsurr] = ll(parameters,dataToFit,PriorMean,PriorInverseCovariance,doPrior,otherOptions)
% 
% where otherOptions.generatesurrogatedata is a binary flag defining whether to apply the prior, and
% doGenerate is a flag defining whether surrogate data (output in asurr) is
% generated. 
% 
% name: names of model likelihood function in folder models
% npar: number of paramters for each 
% parnames: names of parameters for plotting
% partransform: what to do to the (transformed) parameter estimates to transform them into the
% parameters
%  
% SPECIFIC INFO: 
% 
% Simple Rescorla-Wagner model wiht learning rate and beta parameter 

% add models folder to path 
modelDir = 'mBasicRescorlaWagner';
addpath(modelDir);

i=0; 

i=i+1; 
model(i).descr = 'RW model';
model(i).name = 'llrw';			
model(i).npar = 2;
model(i).parnames = {'\beta','\alpha'};
model(i).parnames_untr = {'log \beta','siginv \alpha'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))'};


nModls = i; 
fprintf('%i models loaded from folder %s \n',i,modelDir);
