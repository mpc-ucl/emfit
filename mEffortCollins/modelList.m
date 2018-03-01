function model = modelList; 
% 
% A file like this should be contained in each model class folder and list the
% models to be run, together with some descriptive features. 
% 
% GENERAL INFO: 
% 
% list the models to run here. The models must be defined as likelihood
% functions in the models folder. They must have the form: 
% 
%    [l,dl,dsurr] =
%    ll(parameters,dataToFit,PriorMean,PriorInverseCovariance,doPrior,otherOptions)
% 
% where otherOptions.generatesurrogatedata is a binary flag defining whether to
% apply the prior, and doGenerate is a flag defining whether surrogate data
% (output in asurr) is generated. 
% 
% name: names of model likelihood function in folder models npar: number of
% paramters for each parnames: names of parameters for plotting partransform:
% what to do to the (transformed) parameter estimates to transform them into the
% parameters
%  
% SPECIFIC INFO: 
% 
% This contain models for the effort task by Gold, J. M.; Strauss, G. P.; Waltz,
% J. A.; Robinson, B. M.; Brown, J. K. & Frank, M. J. Negative symptoms of
% schizophrenia are associated with abnormal effort-cost computations. Biol
% Psychiatry, 2013, 74, 130-136
% 
% 
% Quentin Huys 2018 qhuys@cantab.net

i=0; 

i=i+1; 
model(i).descr = 'Effort sensitivity only. This model disregards the reward information and provides a constant bias as a function of the effort associated with the choices.';
model(i).name = 'lleffort';			
model(i).npar = 1;
model(i).parnames = {'\theta'};
model(i).parnames_untr = {'\theta'};
model(i).partransform = {'@(x)x'};

i=i+1; 
model(i).descr = 'Linear effort and reward sensitivity. This model contains an effort and reward sensitivity parameter.';
model(i).name = 'llrewardeffort';			
model(i).npar = 2;
model(i).parnames = {'\beta_{rew}','\beta_{effort}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{effort}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)'};


nModls = i; 
fprintf('%i models in model list\n',i);
