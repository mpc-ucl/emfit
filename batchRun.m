% Perform batch EM inference 
% 
% Quentin Huys, 2018 qhuys@cantab.net
% 
%==============================================================

clear all; 

%--------------------------------------------------------------
% 
% MODELCLASS can take on values
% 
% 	mAffectiveGoNogo
% 	mBasicRescorlaWagner
% 	mProbabilisticReward
% 	mTwostep
%  mPruning 

modelClassToFit = 4; 
 
modelClass{1} = 'mBasicRescorlaWagner';
modelClass{2} = 'mAffectiveGoNogo';
modelClass{3} = 'mProbabilisticReward';
modelClass{4} = 'mTwostep';	% plots... 
% modelClass{5} = 'mPruning';
% modelClass{6} = 'mEffort';
cleanpath(modelClass);

%--------------------------------------------------------------
% add folder to path and load the models 
 
addpath('lib'); 
addpath(genpath(modelClass{modelClassToFit}));
models=modelList; 

%--------------------------------------------------------------
% select models to actually fit and run 
whichinf = [2]; %
models = models(whichinf);

%--------------------------------------------------------------
% load data 
% 
% see the dataformat.txt files in the model folders for instructions on how the
% data contained Data should be formatted for fitting. For demo purposes, we can
% also generate some example surrogate data: 

Data=generateExampleDataset(30); 

%--------------------------------------------------------------
% batchModelFit(Data,modelClass,whichinf,options)
options.checkgradients = 1;			% check gradients of models? 
options.bsub = 0; 						% submit to bsub? 
batchModelFit(Data,models,options); 

%--------------------------------------------------------------
% perform model comparison 
% bestmodel = batchModelComparison(Data,models);

%--------------------------------------------------------------
% generate surrogate data 
% batchGenerateSurrogateData(Data,models);

%--------------------------------------------------------------
% plot surrogate data (specific to model)
load fitResults/SurrogateData; 
surrogateDataPlots(Data,models,SurrogateData,bestmodel)

