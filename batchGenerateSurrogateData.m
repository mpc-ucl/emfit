function SurrogateData = batchGenerateSurrogateData(Data,models,opt)
% 
% SurrogateData = batchModelFit(Data,modelsClass,options)
%
% Generate surrogate experimental data from fitted models contained in modelList
% to compare to original data. 
% 
% options.nSamples	how many samples per subject to generate (default: 100)
% 
% Quentin Huys, 2018 qhuys@cantab.net

if ~exist('opt'); opt = struct; end
if ~isfield(opt,'nSamples'); 		opt.nSamples = 100; 					end;
if ~isfield(opt,'resultsDir'); 	opt.resultsDir= [pwd filesep 'resultsDir'];end;

nModls = length(models);
Nsj = length(Data);

for mdl = 1:nModls

	% load EM-MAP parameters 
	try 
		R.(models(mdl).name) = load(sprintf('%s/%s',opt.resultsDir,models(mdl).name));
		par = R.(models(mdl).name).E; 
	catch 
		fprintf('No fits for model %s found, not surrogate data generated\n',models(mdl).name);
		return 
	end

	fstr=str2func(models(mdl).name);	% turn variable into function call 
	doprior = 0; 							% add a prior for regularization or not
	mu = [];									% prior mean 
	nui = []; 								% prior variance 
	llopt.generatesurrogatedata=1;	% whether to generate surrogate data - here not 

	for sj=1:Nsj;
		fprintf('generating data from model %s.m for subject %i\r',models(mdl).name,sj);
		clear dsurr;
		parfor ns=1:opt.nSamples
			[foo,foo,dsurr(ns)] = fstr(par(:,sj),Data(sj),mu,nui,doprior,llopt); 
		end
		SurrogateData(sj).(models(mdl).name) = dsurr; 
	end
	save([opt.resultsDir '/SurrogateData.mat'],'R', 'SurrogateData');
	fprintf('\n')
end

