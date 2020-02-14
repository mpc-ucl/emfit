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
model(i).descr = 'DDM in combination with a constant model. This model contains one parameter determining the drift rate and parameters for boundary, starting point and non-decision time.';
model(i).name = 'llconstantDDM';				
model(i).npar = 4;
model(i).parnames = {'starting point', 'boundary','theta', 'nonDecisionTime'};
model(i).parnames_untr = {'sig starting pont','log boundary', 'theta', 'sig nonDecisionTime'};
model(i).partransform = {'@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)x', '@(x)1./(1+exp(-x))'};
% % 
% % 
% i=i+1;
% model(i).descr = 'DDM in combination with linear effort and reward sensitivity. This model contains an effort and reward sensitivity parameter determining the drift rate and parameters for boundary and non-decision time.';
% model(i).name = 'llreweffscalingDDMB';				
% model(i).npar = 4;
% model(i).parnames = {'boundary','rew','effort', 'nonDecisionTime'};
% model(i).parnames_untr = {'log boundary', 'log rew','log eff', 'sig nonDecisionTime'};
% model(i).partransform = {'@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)', '@(x)1./(1+exp(-x))'};

% 
i=i+1;
model(i).descr = 'DDM in combination with linear effort and reward sensitivity. This model contains an effort and reward sensitivity parameter determining the drift rate and parameters for boundary, starting point and non-decision time.';
model(i).name = 'llreweffscalingDDMBSP';				
model(i).npar = 5;
model(i).parnames = {'starting point', 'boundary','rew','effort', 'nonDecisionTime'};
model(i).parnames_untr = {'sig starting pont','log boundary', 'log rew','log eff', 'sig nonDecisionTime'};
model(i).partransform = {'@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)', '@(x)1./(1+exp(-x))'};
% % % 
% i=i+1;
% model(i).descr = 'DDM in combination with linear effort and reward sensitivity and probability for post-decision wavering pswitch. This model contains an effort and reward sensitivity parameter determining the drift rate and parameters for boundary, starting point, non-decision time and pswitch.';
% model(i).name = 'llreweffscalingDDMBSPPSwitchEmfit';				
% model(i).npar = 6;
% model(i).parnames = {'starting point', 'boundary','rew','effort', 'nonDecisionTime','pswitch'};
% model(i).parnames_untr = {'sig starting pont','log boundary', 'log rew','log eff', 'sig nonDecisionTime','sig pswitch'};
% model(i).partransform = {'@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)', '@(x)1./(1+exp(-x))', '@(x)1./(1+exp(-x))'};
% 
% i=i+1;
%  model(i).descr = 'DDM in combination with linear effort and reward sensitivity. This model contains an effort and reward sensitivity parameter determining the drift rate and parameters for starting boundary, scaling of boundary,  starting point and non-decision time.';
% model(i).name = 'llreweffscalingDDMBScaledSP';				
% model(i).npar = 6;
% model(i).parnames = {'starting point', 'starting boundary','rew','effort', 'nonDecisionTime', 'boundary scaling'};
% model(i).parnames_untr = {'sig starting pont','log boundary', 'log rew','log eff', 'sig nonDecisionTime','sig boundary scaling'};
% model(i).partransform = {'@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)', '@(x)1./(1+exp(-x))', '@(x)1./(1+exp(-x))'};

 
%  i=i+1;
%  model(i).descr = 'DDM in combination with linear effort and reward sensitivity. This model contains an effort and reward sensitivity parameter determining the drift rate and parameters for starting boundary, scaling of boundary,  starting point and non-decision time and pswitch.';
%  model(i).name = 'llreweffscalingDDMBScaledSPPSwitch';				
%  model(i).npar = 7;
%  model(i).parnames = {'starting point', 'starting boundary','rew','effort', 'nonDecisionTime', 'boundary scaling', 'pswitch'};
%  model(i).parnames_untr = {'sig starting pont','log boundary', 'log rew','log eff', 'sig nonDecisionTime','sig boundary scaling','sig pswitch'};
%  model(i).partransform = {'@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)', '@(x)exp(x)', '@(x)1./(1+exp(-x))', '@(x)1./(1+exp(-x))', '@(x)1./(1+exp(-x))'};



nModls = i; 
fprintf('%i models in model list\n',i);
