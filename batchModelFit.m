function batchEMFit(Data,models,opt)
% 
% batchModelFit(Data,models,options)
% 
% Fit all models contained in MODELCLASS to DATA. 
% 
% MODELCLASS can take on values: 
% 
% 'mBasicRescorlaWagner';			% basic Rescorla-Wagner example 
% 'mAffectiveGoNogo';				% Guitart et al. 2012 
% 'mProbabilisticReward';			% Huys et al., 2013 
% 'mTwostep';							% Daw et al., 2011 
% 'mEffortCollins';	 				% Gold et al., 2013 
% 'mPruning'; 						   % Lally et al., 2017 
% 
% The optional parameter options.CHECKGRADIENTS can be used to check the gradients of
% the provided likelihood files, and options.maxit to limit the number of EM
% iterations. 
%  
% Quentin Huys, 2018 qhuys@cantab.net

% set options 
if ~exist('opt'); opt = struct; end
if ~isfield(opt,'checkgradients'); 	opt.checkgradients = 0; end
if ~isfield(opt,'maxit');           opt.maxit = [];         end

%--------------------------------------------------------------
% Check all the gradients 
%--------------------------------------------------------------
if opt.checkgradients
	for mdl = 1:length(models)
		fprintf('Checking gradient for model %s.m\n',models(mdl).name);
		regressors = cell(models(mdl).npar,1); 		% set up empty regressor cell structure 
		emfit(models(mdl).name,Data,models(mdl).npar,[],[],1); 
		pause 
	end
end

%--------------------------------------------------------------
% EM inference 
%--------------------------------------------------------------

for mdl = 1:length(models)
	fprintf('Performing EM fit %s.m model \n',models(mdl).name);
	regressors = cell(models(mdl).npar,1); 		% set up empty regressor cell structure 
	savestr = sprintf('%s/%s',opt.resultsDir,models(mdl).name);
	[params,var,alpha,stats,bf,fitparams] = emfit(models(mdl).name,Data,models(mdl).npar,[],[],[],[],opt.maxit,[],savestr); 
		
end


