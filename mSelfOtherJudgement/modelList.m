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
% This contain models for the self/other judgement task. 
% 
% Quentin Huys 2021

i=0; 

i=i+1; 
model(i).descr = 'simple model with self and other bias towards accepting positive words';
model(i).name = 'llb';			
model(i).npar = 2;
model(i).parnames = {'bias self','bias other'};
model(i).parnames_untr = {'b_s','b_o'};
model(i).partransform = {'@(x)x','@(x)x'};

i=i+1; 
model(i).descr = 'RW model, restart for each avatar and two self/other bias parameters';
model(i).name = 'll6avb';			
model(i).npar = 4;
model(i).parnames = {'\rho','\alpha','bias self','bias other'};
model(i).parnames_untr = {'log \rho','siginv \alpha','b_s','b_o'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with double update, restart for each avatar and two self/other bias parameters';
model(i).name = 'lld6avb';			
model(i).npar = 4;
model(i).parnames = {'\rho','\alpha','bias self','bias other'};
model(i).parnames_untr = {'log \rho','siginv \alpha','b_s','b_o'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x'};

nModls = i; 
fprintf('%i models in model list\n',i);
