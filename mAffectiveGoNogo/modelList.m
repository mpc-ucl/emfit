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
% This contain models for the affective Go/Nogo task as in: 
% 
% Guitart-Masip M, Huys QJM, Fuentemilla L, Dayan P, Duezel E and Dolan RJ
% (2012): Go and nogo learning in reward and punishment: Interactions between
% affect and effect.  Neuroimage 62(1):154-66 
% 
% Over time, many more variations than models reported in this paper have been
% built and these have been included here. The original set of models were: 
% 
% llba
% llbax
% llbaxq
% llbaxb
% ll2baxb
% llbaepxb
%
% Quentin Huys 2018 qhuys@cantab.net

i=0; 

i=i+1; 
model(i).descr = 'RW model';
model(i).name = 'llba';			
model(i).npar = 2;
model(i).parnames = {'\beta','\alpha'};
model(i).parnames_untr = {'log \beta','siginv \alpha'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise';
model(i).name = 'llbax';			
model(i).npar = 3;
model(i).parnames = {'\beta','\alpha','\gamma'};
model(i).parnames_untr = {'log \beta','siginv \alpha','siginv \gamma'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'RW model with constant bias towards one action ';
model(i).name = 'llbab';			
model(i).npar = 3;
model(i).parnames = {'\beta','\alpha','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','b'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise and with constant bias. ';
model(i).name = 'llbaxb';			
model(i).npar = 4;
model(i).parnames = {'\beta','\alpha','\gamma','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','siginv \gamma','b'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise and with initial bias. ';
model(i).name = 'llbaqx';			
model(i).npar = 4;
model(i).parnames = {'\beta','\alpha','q0','\gamma'};
model(i).parnames_untr = {'log \beta','siginv \alpha','q0','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise and separate reward and loss sensitivities';
model(i).name = 'll2bax';			
model(i).npar = 4;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha','\gamma'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha','siginv \gamma'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias and separate reward and loss sensitivities';
model(i).name = 'll2baxb';			
model(i).npar = 5;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha','\gamma','b'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias and separate reward and loss sensitivities';
model(i).name = 'll2baxb2q';			
model(i).npar = 7;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha','\gamma','b','q0win','q0loss'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha','siginv \gamma','bias','q0_{win}','q0_{loss}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise, initial bias and separate reward and loss sensitivities';
model(i).name = 'll2baqx';			
model(i).npar = 5;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha','q0','\gamma'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha','q0','siginv \gamma'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'RW model with constant bias towards one action, irreducible noise and positive Pavlovian bias parameter';
model(i).name = 'llbaepxb';			
model(i).npar = 5;
model(i).parnames = {'\beta','\alpha','\pi','\gamma','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi','siginv \gamma','b'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with initial bias towards one action, irreducible noise and positive Pavlovian bias parameter';
model(i).name = 'llbaepqx';			
model(i).npar = 5;
model(i).parnames = {'\beta','\alpha','\pi','q0','\gamma'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi','q0','siginv \gamma'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)x','@(x)1./(1+exp(-x))'};

i=i+1; 
model(i).descr = 'RW model with constant bias towards one action, irreducible noise, positive Pavlovian bias parameter and tow different q0 values for win and loss stimuli ';
model(i).name = 'llbaepxb2q';			
model(i).npar = 7;
model(i).parnames = {'\beta','\alpha','\pi','\gamma','b','q0win','q0loss'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi','siginv \gamma','b','q0_{win}','q0_{loss}'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with constant bias towards one action, irreducible noise, separate positive Pavlovian bias parameter for rewards and losses, and two different q0 values for win and loss stimuli ';
model(i).name = 'llba2epxb2q';			
model(i).npar = 8;
model(i).parnames = {'\beta','\alpha','\pi_{rew}','\pi_{loss}','\gamma','b','q0_{win}','q0_{loss}'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi_{win}','log \pi_{loss}','siginv \gamma','b','q0_{win}','q0_{loss}'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, joint reward/loss sensitivity, and separate positively constrained Pavlovian bias parameters for rewards and losses. ';
model(i).name = 'llba2epxb';			
model(i).npar = 6;
model(i).parnames = {'\beta','\alpha','\pi_{rew}','\pi_{loss}','\gamma','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, joint reward/loss sensitivity, and separate positively constrained Pavlovian bias parameters for conflict and non-conflict cues. ';
model(i).name = 'llba2epcxb';			
model(i).npar = 6;
model(i).parnames = {'\beta','\alpha','\pi_{non-confilct}','\pi_{conflict}','\gamma','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi_{non-conflict}','log \pi_{conflict}','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise, constanb bias, separate reward and loss sensitivities, and positive Pavlovian bias parameter. ';
model(i).name = 'll2baepxb';			
model(i).npar = 6;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha','\pi','\gamma','b'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha','log \pi','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, separate reward and loss sensitivities, and separate positively constrained Pavlovian bias parameters for rewards and losses. ';
model(i).name = 'll2ba2epxb';			
model(i).npar = 7;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha','\pi_{rew}','\pi_{loss}','\gamma','b'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

nModls = i; 
fprintf('%i models in model list\n',i);
