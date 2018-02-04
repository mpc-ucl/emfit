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
%    [l,asurr] = ll(parameters,dataToFit,PriorMean,PriorInverseCovariance,doPrior,otherOptions)
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
% This contain models for the two-step Task (Daw et al., 2011 Neuron). 
%  
% Over time, many more variations than models reported in this paper have been
% built and these have been included here. The original model is in llm2b2alr.m
% and contained 7 parameters. 
%  
% Quentin Huys 2018 qhuys@cantab.net

i=0; 

i=i+1; 
model(i).descr = 'Original Daw mixture of model-based and model-free SARSA(lambda) model with sticky repetition factor and separate learning rates for the two levels';
model(i).name = 'llm2b2alr';			
model(i).npar = 7;
model(i).parnames = {'\beta_{1}','\beta_{2}','\alpha_{1}','\alpha_{2}','\lambda','\omega','r'};
model(i).parnames_untr = {'log \beta_1','log \beta_2','siginv \alpha_1','siginv \alpha_2','siginv \lambda','siginv \omega','r'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Pure model-free SARSA(lambda) model with sticky repetition factor and single learning rate';
model(i).name = 'llbmfalr';			
model(i).npar = 4;
model(i).parnames = {'\beta_{MF}','\alpha','\lambda','r'};
model(i).parnames_untr = {'log \beta','siginv \alpha','siginv \lambda','r'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Pure model-based tree-search model with sticky repetition factor and simple RW learning at stage 2.';
model(i).name = 'llbmbar';			
model(i).npar = 4;
model(i).parnames = {'\beta_{MB}','\beta_{MF}','\alpha','r'};
model(i).parnames_untr = {'log \beta_{MB}','log \beta_{MF}','siginv \alpha','r'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Reparametrized mixture of model-based and model-free SARSA(lambda) model with sticky repetition factor and separate learning rates for the two levels. Here, MF and MB values have separate betas, but there is no weight w';
model(i).name = 'll2bmfbmb2alr';			
model(i).npar = 7;
model(i).parnames = {'\beta_{MB}','\beta_{MF1}','\beta_{MF2}','\alpha_{1}','\alpha_{2}','\lambda','r'};
model(i).parnames_untr = {'log \beta_{MB}','log \beta_{MF1}','log \beta_{MF2}','siginv \alpha_1','siginv \alpha_2','siginv \lambda','r'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Reparametrized mixture of model-based and model-free SARSA(lambda) model with sticky repetition factor and separate learning rates for the two levels. Here, MF and MB values have separate betas, but there is no weight w. In addition, MF betas and learning rates at level 1 and 2 are assumed the same.';
model(i).name = 'llbmfbmbalr';			
model(i).npar = 5;
model(i).parnames = {'\beta_{MB}','\beta_{MF}','\alpha','\lambda','r'};
model(i).parnames_untr = {'log \beta_{MB}','log \beta_{MF}','siginv \alpha','siginv \lambda','r'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Pure model-free SARSA(lambda) model with sticky repetition factor and separate learning rates for the two levels';
model(i).name = 'll2bmfalr';			
model(i).npar = 5;
model(i).parnames = {'\beta_{MF1}','\beta_{MF2}','\alpha','\lambda','r'};
model(i).parnames_untr = {'log \beta','siginv \alpha','siginv \lambda','r'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Pure model-free SARSA(lambda) model with sticky repetition factor, separate betas and learning rates for the two levels';
model(i).name = 'll2bmf2alr';			
model(i).npar = 6;
model(i).parnames = {'\beta_{MF1}','\beta_{MF2}','\alpha_{1}','\alpha_{2}','\lambda','r'};
model(i).parnames_untr = {'log \beta_{MF1}','log \beta_{MF2}','siginv \alpha_1','siginv \alpha_2','siginv \lambda','r'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Reparametrized mixture of model-based and model-free SARSA(lambda) model with sticky repetition factor and separate learning rates for the two levels. Here, MF and MB values have separate betas, but there is no weight w. This model additionally contains counterfactional updates in the MF system at level 2 only. ';
model(i).name = 'llbmfbmb2alcr';			
model(i).npar = 7;
model(i).parnames = {'\beta_{MB}','\beta_{MF1}','\beta_{MF2}','\alpha_{1}','\alpha_{2}^{cfct}','\lambda','r'};
model(i).parnames_untr = {'log \beta_{MB}','log \beta_{MF1}','log \beta_{MF2}','siginv \alpha_1','siginv \alpha_2^{cfct}','siginv \lambda','r'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};


nModls = i; 
fprintf('%i models in model list\n',i);
