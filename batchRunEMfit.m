function batchRunEMfit(modelClassToFit,Data,resultsDir,varargin)
% 
% batchRunEMfit(modelClassToFit,pathToData,resultsDir,'key1','val1','key2','val2');
% 
% Performs batch EM inference by first fitting a set of models from each model
% class, plotting the inferred parameters, performing iBIC model comparison,
% generating surrogate data and performing some visual comparisons of the true
% and surrogate data. 
%
% MODELCLASSTOFIT determines which model sets are fitted: 
% 
% 'mBasicRescorlaWagner';			% basic Rescorla-Wagner example 
% 'mAffectiveGoNogo';				% Guitart et al. 2012 
% 'mProbabilisticReward';			% Huys et al., 2013 
% 'mTwostep';							% Daw et al., 2011 
% 'mEffortCollins';	 				% Gold et al., 2013 
% 'mPruning'; 						   % Lally et al., 2017 
% 'mEffortDDM';                     % Berwian et al., 2020
%
% DATA (optional) contains the data.  See the dataformat.txt files in the model
% folders for instructions on how the data contained in DATA should be formatted
% for fitting. For demo purposes, a correct dataset is generated if no data is
% provided. 
% 
% RESULTSDIR (optional) is a path to a directory containing the results. If it
% is not provided then fitResults in the current working directory is used. 
% 
% In addition, the following key-value pairs are accepted: 
% 
% MODELSTOFIT (optional) is a vector with indices of specific models to fit, as
% per the modelList definition in each model folder mXXXX. 
% 
% CHECKGRADIENTS (optional) can be set to 1 to check the gradients of the
% likelihood functions using finite differences. 
% 
% MAXIT (optional) limits the maximal number of EM iterations. 
% 
% Quentin Huys, 2018 qhuys@cantab.net
% 
%==============================================================================

if exist('resultsDir')~=1 || isempty(resultsDir)
	resultsDir = [pwd '/fitResults'];
end
warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir(resultsDir);

% check optional arguments 
validvarargins = {'modelsToFit','checkgradients','maxit','bsub'};
if exist('varargin')
	for k=1:length(validvarargins)
		i = find(cellfun(@(x)strcmpi(x,validvarargins{k}),varargin),1);
		if ~isempty(i); eval([validvarargins{k} '= varargin{i+1};']);end
	end
end

%------------------------------------------------------------------------------
% MODEL CLASS - define which type of model to fit 
 
modelClass{1} = 'mBasicRescorlaWagner';			% basic example 
modelClass{2} = 'mAffectiveGoNogo';					% Guitart et al. 2012 
modelClass{3} = 'mProbabilisticReward';			% Huys et al., 2013 
modelClass{4} = 'mTwostep';							% Daw et al., 2011 
modelClass{5} = 'mEffortCollins';	 				% Gold et al., 2013 
modelClass{6} = 'mPruning'; 							% Lally et al., 2017 
modelClass{7} = 'mEffortDDM'; 							% Berwian et al., 2020 

modelClassToFit = find(cellfun(@(x)strcmp(x,modelClassToFit),modelClass)); 
if isempty(modelClassToFit); error('Model class not found');end

%==============================================================================
% everything below here should just run without alterations 

%------------------------------------------------------------------------------
% get model descriptions 
emfitpath = fileparts(which('batchRunEMfit'));
addpath([emfitpath '/lib']); 
cleanpath(modelClass);									% clean all model paths 
addpath(genpath([emfitpath '/' modelClass{modelClassToFit}]));	% add chosen model path
models=modelList; 										% get complete model list 
if exist('modelsToFit')
	models = models(modelsToFit); 					% select specific models to fit 
end

%------------------------------------------------------------------------------
% generate surrogate data if no data was provided
if ~exist('Data', 'var') || isempty(Data)
	fprintf('No data provided so generating example dataset\n');
	Data=generateExampleDataset(30,resultsDir); 			
end

%------------------------------------------------------------------------------
% fit models using emfit.m
if exist('checkgradients')
	options.checkgradients = checkgradients;		% check gradients of models? 
end
if exist('maxit')
	options.maxit = maxit; 								% limit EM iterations
end
options.resultsDir = resultsDir; 					% directory with results
batchModelFit(Data,models,options); 				% perform the fitting

%------------------------------------------------------------------------------
% perform model comparison 
bestmodel = batchModelComparison(Data,models,resultsDir);

%------------------------------------------------------------------------------
% plot parameters of best model 
batchParameterPlots(Data,models,resultsDir,bestmodel); 

%------------------------------------------------------------------------------
% generate surrogate data 
options.nSamples=100;
batchGenerateSurrogateData(Data,models,options);

%------------------------------------------------------------------------------
% plot surrogate data, compare with real data (specific to model)
load([resultsDir '/SurrogateData']); 
surrogateDataPlots(Data,models,SurrogateData,bestmodel,resultsDir)

fprintf('\n\nAll done!\n\n')
