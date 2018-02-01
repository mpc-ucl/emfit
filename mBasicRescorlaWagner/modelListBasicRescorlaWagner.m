% Model calls 
% 
% list the models to run here. The models must be defined as likelihood functions in
% the models folder. They must have the form: 
% 
%    [l,asurr] = ll(parameters,dataToFit,PriorMean,PriorInverseCovariance,doPrior,otherOptions)
% 
% where doPrior is a binary flag defining whether to apply the prior, and
% otherOptions.generatesurrogatedata is a flag defining whether surrogate data (output in asurr) is
% generated. 
% 
% For each model below, provide the following fiels: 
% 
% model(i).descr: 			description of the model - a string. 
% model(i).name: 				name of model likelihood function 
% model(i).npar: 				number of paramters for each 
% model(i).parnames: 		names of parameters for plotting
% model(i).parnames_untr: 	names of parameters for plotting - with transform 
% model(i).partransform: 	functions to transform inferred values into parameters 

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
