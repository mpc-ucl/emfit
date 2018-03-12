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
% This contain models for the pruning task
% 
%  
% Huys QJM*, Eshel N*, O'Lions E, Sheridan L, Dayan P and Roiser JP (2012):
% Bonsai trees in your head: How the Pavlovian system sculpts goal-directed
% choices by pruning decision trees PLoS Comp Biol 8(3): e1002410 
% 
% Quentin Huys 2018 qhuys@cantab.net


i=0; 

i=i+1; 
model(i).descr = 'Lookahead model. Computes full decision-tree. ';
model(i).name = 'llsb';			
model(i).npar = 1;
model(i).parnames = {'\beta'};
model(i).parnames_untr = {'log \beta'};
model(i).partransform = {'@(x)exp(x)'};

i=i+1; 
model(i).descr = 'Discount model. Flatly discounted decision-tree. ';
model(i).name = 'llsp';			
model(i).npar = 2;
model(i).parnames = {'\beta','\gamma'};
model(i).parnames_untr = {'log \beta','siginv \gamma'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'Pruning model. Discounted decision-tree separately after large losses and other outcomes. ';
model(i).name = 'lls2p';			
model(i).npar = 3;
model(i).parnames = {'\beta','\gamma_G','\gamma_S'};
model(i).parnames_untr = {'log \beta','siginv \gamma_G','siginv \gamma_S'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'Full lookahead but fitting separate reward sensitivities to separate outcomes. ';
model(i).name = 'llsrho';			
model(i).npar = 4;
model(i).parnames = {'\rho_{-70}','\rho_{-20}','\rho_{20}','\rho_{140}'};
model(i).parnames_untr = {'log \rho_{-70}','log \rho_{-20}','log \rho_{20}','log \rho_{140}'};
model(i).partransform = {'@(x)x','@(x)x','@(x)x','@(x)x'};

i=i+1; 
model(i).descr = 'Discount model but fitting separate reward sensitivities to separate outcomes. ';
model(i).name = 'llsrhop';			
model(i).npar = 5;
model(i).parnames = {'\rho_{-70}','\rho_{-20}','\rho_{20}','\rho_{140}','\gamma'};
model(i).parnames_untr = {'\rho_{-70}','\rho_{-20}','\rho_{20}','\rho_{140}','siginv \gamma'};
model(i).partransform = {'@(x)x','@(x)x','@(x)x','@(x)x','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'Pruning model discounting decision-tree separately after large losses and other outcomes, and fitting separate reward sensitivities to separate outcomes. ';
model(i).name = 'llsrho2p';			
model(i).npar = 6;
model(i).parnames = {'\rho_{-70}','\rho_{-20}','\rho_{20}','\rho_{140}','\gamma_G','\gamma_S'};
model(i).parnames_untr = {'\rho_{-70}','\rho_{-20}','\rho_{20}','\rho_{140}','siginv \gamma_G','siginv \gamma_S'};
model(i).partransform = {'@(x)x','@(x)x','@(x)x','@(x)x','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))'};

nModls = i; 
fprintf('%i models in model list\n',i);
