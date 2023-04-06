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
% This contain models for the probabilistic reward task 
% 
% 
% 
% Over time, many more variations than models reported in this paper have been
% built and these have been included here. The original set of models were: 
% 
% llbgeq0.m
% llageqa.m
% ll2bgeq0.m
% llbgelq0.m

%
% Quentin Huys 2018 qhuys@cantab.net

i=0; 

i=i+1; 
model(i).descr = 'Basic stimulus-action Rescorla-Wagner model. This assumes individuals correctly assign the rewards to particular stimulus-action combinations. It has four parameters: a reward learning rate (\alpha) and sensitivity (\beta), an instruction sensitivity (\gamma) and an initial action bias q_0.';
model(i).name = 'llbgeq0';			
model(i).npar = 4;
model(i).parnames = {'\beta','\gamma','\alpha','q_0'};
model(i).parnames_untr = {'log \beta','log \gamma','siginv \alpha','q_0'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Stimulus-action Rescorla-Wagner model with separate sensitivities for reward and non-reward events. This allows for the possibility taht subjects treat non-rewards as actual punishments. This assumes individuals correctly assign the rewards to particular stimulus-action combinations. It has five parameters: a reward learning rate (\alpha), two reward sensitivitites (\beta_reward and \beta_punishment), an instruction sensitivity (\gamma) and an initial action bias q_0.';
model(i).name = 'll2bgeq0';			
model(i).npar = 5;
model(i).parnames = {'\beta_{rew}','\beta_{pun}','\gamma','\alpha','q_0'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{pun}','log \gamma','siginv \alpha','q_0'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Action only model. This assumes individuals only learn about the value of each action, independent of the stimuli. It has four parameters: a reward learning rate (\alpha) and sensitivity (\beta), an instruction sensitivity (\gamma) and an initial action bias q_0. Note that although this model superficially has the same parameters as the model llbgeq0, the meaning of these parameters is different.';
model(i).name = 'llageqa';			
model(i).npar = 4;
model(i).parnames = {'\beta','\gamma','\alpha','q_a'};
model(i).parnames_untr = {'log \beta','log \gamma','siginv \alpha','q_a'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Belief model. As subjects are unsure about the presented stimulus, the might assign rewards to both stimuli, with only a certain preference for the actually presented stimulus. The model has five parameters: a reward learning rate (\alpha), a reward sensitivity \beta, a belief bl, an instruction sensitivity \gamma and an initial action bias q_0.';
model(i).name = 'llbgelq0';			
model(i).npar = 5;
model(i).parnames = {'\beta','\gamma','\alpha','bl','q_0'};
model(i).parnames_untr = {'log \beta','log \gamma','siginv \alpha','siginv bl','q_0'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'Belief model with counterfactual updates. As subjects are unsure about the presented stimulus, the might assign rewards to both stimuli, with only a certain preference for the actually presented stimulus. It additionally performs counterfactual updates. The model has five parameters: a reward learning rate (\alpha), a reward sensitivity \beta, a belief bl, an instruction sensitivity \gamma and an initial action bias q_0.';
model(i).name = 'lldbgelq0';			
model(i).npar = 5;
model(i).parnames = {'\beta','\gamma','\alpha','bl','q_0'};
model(i).parnames_untr = {'log \beta','log \gamma','siginv \alpha','siginv bl','q_0'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};


nModls = i; 
fprintf('%i models in model list\n',i);
