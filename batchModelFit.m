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
	if ~opt.bsub
		[params,var,alpha,stats,bf,fitparams] = emfit(models(mdl).name,Data,models(mdl).npar,[],[],[],[],[],[],savestr); 
	else
		
		% bsub -W 8:00 -o ${outdir}/${i}pre -J ${i}pre ${matlab} -r "preprocPruning($i)"

		% numproc=length(models)
		% 
		% # where to write output files (what matlab would display)
		% batchname='batch' # what the name of the output files will be
		% outdir='/cluster/scratch/huysq/aida/pruningfmri/data/'${batchname}'/'
		% matlab='/cluster/apps/matlab/9.1/x86_64/bin/./matlab -nodisplay -nosplash -singleCompThread'
		% 
		% rm -r $outdir
		% mkdir -p $outdir
		% 
		% for ((i=1;i<=$numproc;i++)); do
		%    # first run the preprocessing
		%    # then run firstlevel analysis after preprocessing is done
		%    # bsub -w "done(${i}pre)" -W 1:00 -o ${outdir}/${i}fir -J ${i}fir ${matlab} -r "firstlevel($i)"
		%    bsub                     -W 1:00 -o ${outdir}/${i}fir -J ${i}fir ${matlab} -r "firstlevel($i)"
		%    waitlist[$i]="-w done(${i}fir)"
		% done
		% 
		% bsub ${waitlist[*]} -W 4:00 -o ${outdir}/secondlevel -J secondlevel ${matlab} -r "secondlevel"

	end
end

% 
% 
