% Perform batch EM inference 
% 
% Quentin Huys, 2018 qhuys@cantab.net
% 
%==============================================================

clear all; 

%--------------------------------------------------------------
% load data 
% 
% see the dataformat.txt files in the model folders for instructions on how the
% data contained Data should be formatted for fitting. 

load data/exampleAffectiveGoNogoData; 

%--------------------------------------------------------------
% 
% MODELCLASS can take on values
% 
% 	mAffectiveGoNogo
% 	mBasicRescorlaWagner
% 	mProbabilisticReward
% 	mTwostep
 
modelClass = 'mAffectiveGoNogo';
modelClass = 'mBasicRescorlaWagner';

%--------------------------------------------------------------
% add folder to path and load the models 
 
% cd /Users/qhuys/git/emfittoolbox/	% cd to emfit folder 
addpath('lib'); 
addpath(genpath(modelClass));
models=modelList; 

%--------------------------------------------------------------
% select models to actually fit and run 
% whichinf = [1 4 7]; %
% models = models(whichinf);

%--------------------------------------------------------------
% batchModelFit(Data,modelClass,whichinf,checkgradients)
batchModelFit(Data,models); 

%--------------------------------------------------------------
% perform model comparison 
bestmodel = batchModelComparison(Data,models);

%--------------------------------------------------------------
% generate surrogate data 
batchGenerateSurrogateData(Data,models);

%--------------------------------------------------------------
% plot surrogate data (specific to model)
load fitResults/SurrogateData; 
surrogateDataPlots(Data,models,SurrogateData,bestmodel)

