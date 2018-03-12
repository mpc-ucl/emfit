%==============================================================================
% 
% Perform batch EM inference 
% 
% Quentin Huys, 2018 qhuys@cantab.net
% 
%==============================================================================

clear all; 

%------------------------------------------------------------------------------
% MODEL CLASS - define which type of model to fit 
 
modelClassToFit = 6; 	% choose one of the classes below 
 
modelClass{1} = 'mBasicRescorlaWagner';			% basic example 
modelClass{2} = 'mAffectiveGoNogo';					% Guitart et al. 2012 
modelClass{3} = 'mProbabilisticReward';			% Huys et al., 2013 
modelClass{4} = 'mTwostep';							% Daw et al., 2011 
modelClass{5} = 'mEffortCollins';	 				% Gold et al., 2013 
modelClass{6} = 'mPruning'; 							% Lally et al., 2017 

%------------------------------------------------------------------------------
% LOAD DATA 
% 
% see the dataformat.txt files in the model folders for instructions on how the
% data contained Data should be formatted for fitting. For demo purposes, we can
% also generate some example surrogate data: 

generateExampleDataset(30); 			% only if you don't have your own data! 
load Data; 


%==============================================================================
% everything below here should just run without alterations 

%------------------------------------------------------------------------------
% get model descriptions 
addpath('lib'); 
cleanpath(modelClass);									% clean all model paths 
addpath(genpath(modelClass{modelClassToFit}));	% add chosen model path
models=modelList; 										% get complete model list 
% models = models([1,3,7]); 							% select specific models to fit 

%------------------------------------------------------------------------------
% fit models using emfit.m
options.checkgradients = 0;							% check gradients of models? 
options.bsub = 0; 										% submit to bsub? - in progress 
batchModelFit(Data,models,options); 				% perform the fitting
batchParameterPlots(Data,models,bestmodel);		% plot parameters 

%------------------------------------------------------------------------------
% perform model comparison 
bestmodel = batchModelComparison(Data,models);

%------------------------------------------------------------------------------
% generate surrogate data 
options.nSamples=100;
batchGenerateSurrogateData(Data,models,options);

%------------------------------------------------------------------------------
% plot surrogate data, compare with real data (specific to model)
load fitResults/SurrogateData; 
surrogateDataPlots(Data,models,SurrogateData,bestmodel)

fprintf('\n\nAll done!\n\n')
