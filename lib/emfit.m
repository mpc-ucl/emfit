function [E,V,alpha,stats,bf,fitparams] = emfit(llfunc,D,Np,varargin); 
%
% [E,V,alpha,stats,bf,fitparams] = EMFIT(llfunc,D,Np,[reg],[Nsample],[docheckgrad],[nograd],[maxit],[dofull],[savestr],[loadstr]); 
%  
% Perform a random-effects fit using expectation-maximimization. 
%
% NOTE: This is in development. The error bars around the group mean are only
% correct for small models with few parameters, not for larger ones. 
% 
% LLFUNC	- this is a string that points to a likelihood function of choices (a
% model). The function must have the form: 
%
%    [l,dl] = llfunc(x,D,musj,nui,doprior,llopt);
% 
% where x are the parameters, D(sj) is the data for subject sj. MUSJ(:,sj) is
% the prior mean parameter for subject sj. If no regressors (see below) are
% included, then this should be the same for all subjects. If a regressor is
% included for a GLM analysis, then each subject's value will depend on this.
% DOPRIOR determines whether a Gaussian prior with mean musj and diagonal
% covariance matrix NU is added to the log likelihood of the choices.  The
% function must output the sum of the total log likelihood of the choice (l),
% and the gradient wrt all the paramters (dl). For example llfunc see llrw.m
% included at the end of this function 
% 
% D contains all the data, with D(sj) containing all the data for subject
% sj=1...Nsj. If D(sj).Nch are the number of observations for subjects sj, then
% bf can return integrated BIC values and integrated Laplacian values (see
% below). 
% 
% NP is th enumber of parameters in the function llfunc. 
%  
% REG is a cell structure that must be the length NP.  For each parameter, REG
% can be defined. For instance, for a model with two parameters, one might want
% to ask whether the second parameter is related to some psychometric variable
% psi at the group level, taking into account how noisy each subject is. To do
% this, define 
%  
%     reg{1} = []; 
%     reg{2} = psi;   
% 
% with psi(i) being the psychometric score for subject i. 
% 
% PARALLEL: If a matlabpool has been opened prior to calling the function it
% will use it (using PARFOR) and hence run faster. Ideally, you should use the
% same number of processors as you have separate subjects in your dataset, or at
% least around 80%. 
%   
% The OUTPUT is as follows: 
%     E        is a matrix of size NpxNsj containing the MAP-EM parameters
%     V        is a matrix of size NpxNsj containing the inverse Hessians around 
%              individual subject parameter estimates 
%     alpha    contains the estimated coefficients (both means and regression 
%              coefficients if REG has been included)
%     stats    contains further stats for each estimated prior parameter, in particular
%              standard error estimates (from finite differences), t and p
%              values, and ML estimates as stats.EML. stats.EX is 1 if it
%              converged; 0 if it reached the MAXIT limit, -2 if a saddle point
%              (with some positive group Hessian diagonals) rather than a maximum was
%              reached, and -1 otherwise. 
%     bf       contains estimates of the quantities necessary to compute Bayes
%    			   factors for model comparison. The integrated, group-level BIC
%    			   bf.ibic simply counts the parameters at the group level, while
%    			   bf.ilap uses the finite difference Hessian around the prior mean
%    			   for a Laplacian estimate, which is typically less conservative. 
%     fitparam If asked for, this outputs various parameters used by emfit for
%              fitting - including the data (this might take up a lot of space!)
%  
% Additional variables can be defined: 
%     
%     NSAMPLES is the number of samples used for integration (default: 2000). 
%     DOCHECKGRAD=1 sets a flag to check the gradients dl provided by llfunc
%     (default 0). NOGRAD has to be set to 1 if no gradients are provided by
%     llfunc. MAXIT is the maximum number of EM iterations (default 500). If this
%     is set to 0, then ML and only a single EM iteration will be computed. The
%     ML parameters are always in stats.EML. If DOFULL is set to zero, only the
%     diagonal of the Hessians are used for EM. This effectively imposes a prior
%     that sets off-diagonal elements to zero, but is no longer recommended.
%     If SAVESTR is provided, intermediate results are saved in this file. 
% 
% Copyright Quentin Huys and Daniel Schad, 2015 
% www.quentinhuys.com/code.html 
% www.quentinhuys.com/pub.html
% qhuys@cantab.net


%=====================================================================================
% setting up 

fprintf('---------------------------------------------------------------------------\n')
fprintf('NOTE: emfit.m is in development. Statistics on MAP and EM-MAP\n');
fprintf('parameters are stable. The error bars around the group mean \n');
fprintf('are only correct for small models. \n')
fprintf('---------------------------------------------------------------------------\n')

fitparams.dx= 0.001; 									% step for finite differences
fitparams.tol = 1e-8; 									% fminunc tolerance 
fitparams.robust = 3; 									% restarts of each individual fit (at least 1)
fprintf('performing %i restarts for every internal fminunc call\n',fitparams.robust);

fstr=str2func(llfunc);									% prepare function string 
Nsj = length(D); sjind=1:Nsj;							% number of subjects 

nargin = length(varargin); 
t=varargin; 
if nargin>0 & ~isempty(t{1}); reg         = t{1}; else reg=cell(Np,1);   end; 
if nargin>1 & ~isempty(t{2}); Nsample     = t{2}; else Nsample = 2000;   end; 
if nargin>2 & ~isempty(t{3}); docheckgrad = t{3}; else docheckgrad = 0 ; end; 
if nargin>3 & ~isempty(t{4}); nograd      = t{4}; else nograd = 0 ;      end; 
if nargin>4 & ~isempty(t{5}); maxit       = t{5}; else maxit = 500;      end; 
if nargin>5 & ~isempty(t{6}); dofull      = t{6}; else dofull = 1;       end; 
if nargin>6 & ~isempty(t{7}); savestr     = t{7}; else savestr = '';     end; 
if nargin>7 & ~isempty(t{8}); loadstr     = t{8}; else loadstr = '';     end; 

% never simulate surrogate data while fitting: 
llopt.generatesurrogatedata=0;

% deal with gradients being provided or not 
if nograd; 													% assume gradients are supplied 
	fminopt=optimoptions('fminunc','display','off','TolFun',fitparams.tol);
	%fminopt = optimoptions(@fminunc,'Display','off','TolFun',fitparams.tol);
else 
	fminopt=optimoptions('fminunc','display','off','GradObj','on','TolFun',fitparams.tol,'algorithm','trust-region');
	%fminopt = optimoptions(@fminunc,'Display','off','GradObj','on','TolFun',fitparams.tol);
	if docheckgrad; 											% if they are, then can check them. 
		fminopt=optimoptions('fminunc','display','off','GradObj','on','DerivativeCheck','on','TolFun',fitparams.tol);
		%fminopt = optimoptions(@fminunc,'Display','off','GradObj','on','DerivativeCheck','on','TolFun',fitparams.tol);
	end
end
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
%warning('off','optim:fminunc:WillRunDiffAlg');

% check if regressors have been provided correctly 
if Np~=length(reg); 
	error('You must provide a regressor cell entry for each parameter.');
elseif ~iscell(reg);
	error('REG must be a cell structure of length Np.');
end

% make regression matrices 
Xreg = repmat(eye(Np),[1 1 Nsj]);
Nreg=0; 
for j=1:Np
	if size(reg{j},2)==Nsj; reg{j} = reg{j}';end
	for k=1:size(reg{j},2)
		Nreg = Nreg+1;
		Xreg(j,Np+Nreg,:) = reg{j}(:,k); 
	end
end
Npall= Np+Nreg; 
coeff_vec = Inf*ones(Npall,1);

% load & continue previous fit or set up things for a new fit 


try 																% Try continuing previous fit 
	eval(['load ' loadstr ' E V alpha emit fitparams musj nui stats ']);
	nu = inv(nui);
	fprintf('Loaded %s, continuing fit will save as %s\n',loadstr,savestr); 
catch 															% else set things up for new fit 
	alpha = zeros(Npall,1); 								% group regression coefficients 
	for sj=1:Nsj
		musj(:,sj) = Xreg(:,:,sj)*alpha; 				% individual subject means
	end
	nui = 0.1*eye(Np); nu = inv(nui);					% prior variance over all params 
	E = zeros(Np,Nsj) + sqrtm(nu)*randn(Np,Nsj);		% random initial individual subject parameter estimates
	emit=0;stats.ex=-1;
	fprintf('No old fit loaded, starting new fit');if ~isempty(savestr); fprintf(', will save as %s\n',savestr); else fprintf('\n');end
end

% check gradients - always with prior 
if docheckgrad; checkgrad(func2str(fstr),randn(Np,1),.001,D(1),musj(:,1),nui,1,llopt), return; end

% prepare for multiple restarts 
sjind_rs = repmat(sjind,[fitparams.robust,1]) + repmat((0:fitparams.robust-1)'*Nsj,[1,Nsj]);
sjind_rs = sjind_rs'; sjind_rs = sjind_rs(:)';

%=====================================================================================
fprintf('\nStarting EM estimation');
	
PLold= -Inf; nextbreak=0;
while 1;emit=emit+1; t0=tic;

	% E step...........................................................................

	% start with ML estimates & set things up for multiple restarts 
	if emit==1; doprior=0; sjind_pf = sjind; else doprior=1; sjind_pf = sjind_rs; end	
	if emit==2; 
		E = zeros(Np,Nsj) + sqrtm(nu)*randn(Np,Nsj);		% random initial individual subject parameter estimates
	end
	% main loop over subjects 
	for sj=sjind_pf; tt0=tic; 
		sk = mod(sj-1,Nsj)+1; 							% current subject
		rs = ceil(sj/Nsj);								% current restart for that subject 
		est=[]; fval=[]; ex=-1; hess=[]; nfc=1; 
		while 1 												% have to deal with poor convergence 
			init = E(:,sk);
			if rs>1 | nfc>1; init = init+ nfc*.1*real(sqrtm(nu))*randn(Np,1); end % add noise for next attempt
			sto = 1; 
            while sto == 1
                try
                    %init(end) = -1.6; 
                    [est(:,nfc),fval(nfc),ex(nfc),foo,foo,hess(:,:,nfc)] = fminunc(@(x)fstr(x,D(sk),musj(:,sk),nui,doprior,llopt),init,fminopt);
                    sto = 0;
                catch
                    F = zeros(Np,Nsj) + sqrtm(nu)*randn(Np,Nsj);
                    E(:,sk) = F(:,sk); 
                    init = E(:,sk);
                    init = init+ nfc*.1*real(sqrtm(nu))*randn(Np,1);
                end
            end
                
			if  any(ex(nfc)==[1:3]) | emit==1 | nfc==3; break;end
			if ~any(ex(nfc)==[1:3]) ; fprintf('sj %i, rep %i convergence failure %i exit status %i\n',sk,rs,nfc,ex); end
			nfc=nfc+1; 
		end
		[foo,best] = min(fval);							% take fit that led to best function value among those that converged
		tE(:,sj)		= est(:,best);						% parameter estimates 
		tW(:,:,sj) 	= pinv(full(hess(:,:,best)));	% covariance matrix around parameter estimate
		tV(:,sj)		= diag(tW(:,:,sj));				% diagonal undertainty around parameter estimate
		tPL(sj) 		= fval(best);						% posterior likelihood 
		ttt(sj) 		= toc(tt0); 						% time info 
		tEx(sj) 		= ex(best); 						% exit codes - for debugging
	end
	% choose best of restarts 
	if emit==1; best = ones(1,Nsj); 					% for first (ML) iteration 
	else;       [PL,best] = min(reshape(tPL',[Nsj fitparams.robust])');	% unwrap multiple restarts
	end
	for sj=1:Nsj; 
		sk = sj+(best(sj)-1)*Nsj;						% take fit that led to best function value over the restarts
		E(:,sj)	= tE(:,  sk);                    % parameter estimates 
		W(:,:,sj)= tW(:,:,sk);                    % covariance matrix around parameter estimate
		V(:,sj)	= tV(:,  sk);                    % diagonal undertainty around parameter estimate
		PL(sj) 	= tPL(   sk);                    % posterior likelihood 
		tt(sj) 	= ttt(   sk);                    % time info 
		Ex(sj) 	= tEx(   sk);                    % exit codes - for debugging
	end
 
	fprintf('\nEMit=%i/%i PL=%.2f loop=%.2gs mean/subj=%.2gs longest=%.2gs parfor speedup=%.2g',emit,maxit,sum(PL),toc(t0),mean(tt),max(tt),mean(tt)*Nsj*fitparams.robust/toc(t0))
	if emit==1; stats.EML = E;   stats.VML = V; end
	if emit==2; stats.EMAP0 = E; stats.VMAP0 = V; end
	if nextbreak==1; break;end

	% M step...........................................................................

	if emit> 1; 							% only update priors after first MAP iteration
		while 1								% iterate over prior mean and covariance updates until convergence
			% prior mean update - note prior mean is different for each subject if
			% GLM is included 
			alpha_old = alpha; 
			ah = zeros(Npall,1); 
			vh = zeros(Npall);
			for  sj=1:Nsj
				ah = ah + Xreg(:,:,sj)'*nui*E(:,sj); 
				vh = vh + Xreg(:,:,sj)'*nui*Xreg(:,:,sj);
			end
			%alpha = pinv(vh)*ah; 
			alpha = vh\ah; 
			for sj=1:Nsj
				musj(:,sj) = Xreg(:,:,sj)*alpha; 
			end

			% prior covariance update 
			if ~dofull 											% use diagonal prior variance 
				nunew = diag(sum(E.^2 + V - 2*E.*musj + musj.^2,2)/(Nsj-1)); 
			else													% use full prior covariance matrix 
				nunew = zeros(Np); 
				for sk=sjind;
					nunew = nunew + E(:,sk)*E(:,sk)' - E(:,sk)*musj(:,sk)'-musj(:,sk)*E(:,sk)' + musj(:,sk)*musj(:,sk)' + W(:,:,sk); 
				end
				nunew = nunew/(Nsj-1); 					
			end
			if det(nunew)<0; fprintf('negative determinant - nufull not updated'); 
			else           ; nu = nunew;
			end
			nui = pinv(nu); 
			if norm(alpha-alpha_old)<1e7; break;end 
		end
	end

	% check for convergence of EM procedure or stop if only want ML / MAP0..............
	if maxit==1 | (maxit==2 & emit==2); break; end		% stop if only want ML or ML& MAP0
	if emit>1;if abs(sum(PL)-PLold)<1e-3;nextbreak=1;stats.ex=1; fprintf('...converged');end;end
	if emit>=maxit; nextbreak=1;stats.ex=0;fprintf('...maximum number of EM iterations reached');end
	PLold=sum(PL);
	stats.diagnostics.sumPL(emit) = sum(PL); 
	stats.diagnostics.mu(:,emit) = mean(musj,2); 
	stats.diagnostics.nu(:,:,emit) = nu; 
	stats.diagnostics.exitcodes(:,emit) = Ex;
	stats.diagnostics.E(:,:,emit) = E; 
	stats.diagnostics.V(:,:,emit) = V; 
	stats.diagnostics.W(:,:,:,emit) = W; 
	if length(savestr)>0; eval(['save ' savestr ' E V alpha stats emit musj nui fitparams']);end
end
stats.PL = PL; 
stats.subjectmeans= musj; 
stats.groupvar= nu; 
stats.alpha = alpha; 
if emit>=3
	stats.E_EMMAP = E; 
	stats.V_EMMAP = V; 
end

%=====================================================================================
fprintf('\nComputing individual subject BIC values');

% Check if observation count has been provided
if ~isfield(D,'Nch');
	error('D.Nch needs to contain the number of observations for each subject');
end

for sk=sjind
	stats.LL(sk)  = fstr(E(:,sk),D(sk),musj(:,sk),nui,0,llopt);
	bf.bic(sk) =  -2*stats.LL(sk)   + Np*log(D(sk).Nch);
end
	
if maxit<=2; return; end								% end here if only want ML or ML & MAP0

%=====================================================================================
fprintf('\n');

dx = fitparams.dx;
oo = ones(1,Nsample);
LLi=zeros(Nsample,1);
Hess = zeros(Npall,Npall,Nsj);
for sk=sjind;
	fprintf('\rSampling subject %i ',sk)

	% estimate p(choices | prior parameters) by integrating over individual parameters
	muo = musj(:,sk)*oo; 
	es = sqrtm(nu)*randn(Np,Nsample)+muo; 
	parfor k=1:Nsample;
		LLi(k) = fstr(es(:,k),D(sk),musj(:,sk),nui,0,llopt); 
	end
	lpk0 = max(-LLi);
	pk0 = exp(-LLi-lpk0);
	bf.iL(sk) = log(mean(exp(-LLi-lpk0)))+lpk0;	% integrated likelihood 

	bf.SampleProbRatio(sk)=(sum(-LLi)/-stats.PL(sk))/Nsample; 
	bf.EffNsample(sk)=sum(exp(-LLi/D(sk).Nch));	% eff # samples

	% shift samples to get gradients and (full) Hessian around prior parameters 
	des = es-muo; 
	err = sum(des.*(nui*des),1); 
	for l=1:Npall
		% shift *samples* by +2*dx 
		foo = alpha; foo(l) = foo(l)+2*dx; mud = Xreg(:,:,sk)*foo; 
		desshift = es - mud*oo; 
		lw = -1/2*sum(desshift.*(nui*desshift),1) - (-1/2*err); 
		w = exp(lw'); w = w/sum(w); 
		lldm = log(pk0'*w)+lpk0;

		% shift *samples* by +dx 
		foo = alpha; foo(l) = foo(l)+dx; mud = Xreg(:,:,sk)*foo; 
		desshift = es - mud*oo; 
		lw = -1/2*sum(desshift.*(nui*desshift),1) - (-1/2*err); 
		w = exp(lw'); w = w/sum(w); 
		lld(l,1) = log(pk0'*w)+lpk0;
		bf.EffNsampleHess(l,l,sk) = 1/sum(w.^2);

		% finite difference diagonal Hessian elements
		Hess(l,l,sk) = (lldm+bf.iL(sk)-2*lld(l,1))/dx^2; 

		% now compute off-diagonal terms
		for ll=1:l-1
			% again shift samples, but along two dimensions 
			foo  = alpha;  foo(l) = foo(l)+dx; foo(ll) = foo(ll)+dx; mud = Xreg(:,:,sk)*foo;
			desshift = es - mud*oo; 
			lw = -1/2*sum(desshift.*(nui*desshift),1) - (-1/2*err); 
			w = exp(lw'); w = w/sum(w); 
			lldd(l,ll) = log(pk0'*w)+lpk0;  % off-diagonal second differential 
			bf.EffNsampleHess(l,ll,sk) = 1/sum(w.^2);	% eff # samples 
			bf.EffNsampleHess(ll,l,sk) = 1/sum(w.^2);

			% off-diagonal Hessian terms 
			Hess(l,ll,sk) = (lldd(l,ll) - lld(l,1) - lld(ll,1) + bf.iL(sk))/dx^2; 
			Hess(ll,l,sk) = Hess(l,ll,sk);
		end
	end
	if any(any(bf.EffNsampleHess(:,:,sk)<50)) 
		warning('Warning: Less than 50 effective samples - dimensionality prob.  too high!');
	end
end
fprintf('...done ')

stats.individualhessians = Hess; 
stats.groupmeancovariance = - pinv(sum(Hess,3))*Nsj/(Nsj-1);
saddlepoint=0;
if any(diag(stats.groupmeancovariance)<0); 
	warning('Negative Hessian, i.e. not at maximum - try running again, increase MAXIT if limit reached')
	stats.ex=-2; 
	saddlepoint=1;
end

% compute t and p values comparing each parameter to zero 
if ~saddlepoint												% can estimate p values 
	stats.groupmeanerr	= sqrt(diag(stats.groupmeancovariance));
	stats.tval				= alpha./stats.groupmeanerr; 
	stats.p					= 2*tcdf(-abs(stats.tval), Nsj-Npall);
else																% can't estimate p values 
	stats.groupmeanerr	= NaN*alpha; 
	stats.tval				= NaN*alpha; 
	stats.p					= NaN*alpha; 
end

%=====================================================================================
fprintf('\nComputing iBIC and iLAP')

Nch=0; for sj=sjind; Nch = Nch + D(sj).Nch;end
bf.ibic =  -2*(sum(bf.iL) - 1/2*(2*Np+Nreg)*log(Nch));
bf.ilap =  -2*(sum(bf.iL) - 1/2*   Np      *log(Nch)  + 1/2*log(det(stats.groupmeancovariance)));

if nargout>=6 
	%=====================================================================================
	fitparams.likelihoodfunction=llfunc; 
	fitparams.reg=reg; 
	fitparams.Nsample=Nsample; 
	fitparams.docheckgrad=docheckgrad; 
	fitparams.nograd=nograd; 
	fitparams.maxit=maxit; 
	fitparams.dofull=dofull; 
	fitparams.D=D; 
	fitparams.Np=Np; 
end

if length(savestr)>0; 
	fprintf('\nSaving data to %s ',savestr)
	eval(['save ' savestr ' E V alpha stats emit musj nui bf fitparams']);
end
	
fprintf('\nDone\n******************\n ')
return 

