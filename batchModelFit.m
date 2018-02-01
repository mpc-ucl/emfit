function batchEMFit(Data,models,opt)
% 
% batchModelFit(Data,models,options)
% 
% Fit all models contained in MODELCLASS to DATA. 
% 
% MODELCLASS can take on values: 
% 
% 	 mAffectiveGoNogo
% 	 mBasicRescorlaWagner
% 	 mProbabilisticReward
% 	 mTwostep
% 
% The optional parameter options.CHECKGRADIENTS can be used to check the gradients of
% the provided likelihood files 
%  
% Quentin Huys, 2018 qhuys@cantab.net

% set options 
if ~exist('opt'); opt = struct; end
if ~isfield(opt,'checkgradients'); 	opt.checkgradients = 0; 				end

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
	savestr = sprintf('fitResults/%s',models(mdl).name);
	[params,var,alpha,stats,bf,fitparams] = emfit(models(mdl).name,Data,models(mdl).npar,[],[],[],[],[],[],savestr); 
end


