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

% NEED TO ADD
% 
% Quentin Huys 2018 qhuys@cantab.net

i=0; 
% 
% i=i+1;
% model(i).descr = 'Trade off of reward and effort';
% model(i).name = 'llreweffch';				
% model(i).npar = 2;
% model(i).parnames = {'reward sensitivity', 'effort sensitivity'};
% model(i).parnames_untr = {'log reward', 'log effort'};
% model(i).partransform = {'@(x)exp(x)','@(x)exp(x)'};

% i=i+1;
% model(i).descr = 'Trade off of reward and effort and increase in effort with number of trials';
% model(i).name = 'llreweffchfat';				
% model(i).npar = 2;
% model(i).parnames = {'reward sensitivity', 'effort sensitivity'};
% model(i).parnames_untr = {'log reward', 'log effort'};
% model(i).partransform = {'@(x)exp(x)','@(x)exp(x)'};
% 
% i=i+1;
% model(i).descr = 'Trade off of reward and effort and rate';
% model(i).name = 'llreweffratech';				
% model(i).npar = 3;
% model(i).parnames = {'reward sensitivity', 'effort sensitivity', 'rate sensitivity'};
% model(i).parnames_untr = {'log reward', 'log effort', 'log rate'};
% model(i).partransform = {'@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)'};
% 
% 
% i=i+1;
% model(i).descr = 'Trade off of reward and effort and average reward rate and max reward';
% model(i).name = 'llreweffavch';				
% model(i).npar = 4;
% model(i).parnames = {'reward sensitivity', 'effort sensitivity', 'rate sensitivity', 'max reward'};
% model(i).parnames_untr = {'log reward', 'log effort', 'log rate', 'log max reward'};
% model(i).partransform = {'@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)', '@(x)exp(x)'};
% % 
% 

% i=i+1;
% model(i).descr = 'Trade off of reward and effort and rate with decrease in rate with number of trials';
% model(i).name = 'llreweffratefatch';				
% model(i).npar = 3;
% model(i).parnames = {'reward sensitivity', 'effort sensitivity', 'rate sensitivity'};
% model(i).parnames_untr = {'log reward', 'log effort', 'log rate'};
% model(i).partransform = {'@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)'};

i=i+1;
model(i).descr = 'Trade off of reward and effort for two groups';
model(i).name = 'llreweffchRel';				
model(i).npar = 4;
model(i).parnames = {'reward sensitivity 1', 'effort sensitivity 1', 'reward sensitivity 2', 'effort sensitivity 2'};
model(i).parnames_untr = {'log reward 1', 'log effort 1', 'log reward 2', 'log effort 2'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)exp(x)','@(x)exp(x)'};


nModls = i; 
fprintf('%i models in model list\n',i);
