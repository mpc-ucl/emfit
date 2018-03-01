%==============================================================================
% 
% Perform batch EM inference 
% 
% Quentin Huys, 2018 qhuys@cantab.net
% 
%==============================================================================

clear all; 

%------------------------------------------------------------------------------
% model class - define which type of model to fit 
 
modelClassToFit = 1; 
 
modelClass{1} = 'mBasicRescorlaWagner';
modelClass{2} = 'mAffectiveGoNogo';
modelClass{3} = 'mProbabilisticReward';
modelClass{4} = 'mTwostep';
modelClass{5} = 'mEffortCollins';	 
modelClass{6} = 'mPruning';	 % old format, doesn't work with batch yet 

%------------------------------------------------------------------------------
% add folder to path and load the models 
 
addpath('lib'); 
cleanpath(modelClass);
addpath(genpath(modelClass{modelClassToFit}));
models=modelList; 

%------------------------------------------------------------------------------
% select models to actually fit and run 
% whichinf = [1:2]; %
% models = models(whichinf);

%------------------------------------------------------------------------------
% load data 
% 
% see the dataformat.txt files in the model folders for instructions on how the
% data contained Data should be formatted for fitting. For demo purposes, we can
% also generate some example surrogate data: 

Data=generateExampleDataset(30); 

%------------------------------------------------------------------------------
% fit models using emfit.m
options.checkgradients = 0;			% check gradients of models? 
options.bsub = 0; 						% submit to bsub? - in progress 
batchModelFit(Data,models,options); 

%------------------------------------------------------------------------------
% perform model comparison 
bestmodel = batchModelComparison(Data,models);

%------------------------------------------------------------------------------
% generate surrogate data 
options.nSamples=100;
batchGenerateSurrogateData(Data,models,options);

%------------------------------------------------------------------------------
% plot surrogate data (specific to model)
load fitResults/SurrogateData; 
surrogateDataPlots(Data,models,SurrogateData,bestmodel)

