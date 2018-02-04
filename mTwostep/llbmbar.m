function [l,dl,dsurr] = llbmbar(x,D,mu,nui,doprior,options);
%
% [l,dl,dsurr] = llbmbar(x,D,mu,nui,doprior,options);
%
% Fit only the model-based component of the two-step task by Daw et al. 2011
% 
% X are the parameters for which to evaluate the likelihood D contains the data
% (see dataformat.txt or generate some data using generateExampleDataset.m).  MU
% is the prior mean and NUI is the prior inverse covariance matrix.  The DOPRIOR
% flag defines whether to apply a prior (1) or not (0). The variable
% OPTIONS.generatesurrogatedata defines whether to return the likelihood of data
% D (0) or whether to generate new surrogate data from the given parameters (1). 
%
% Quentin Huys, 2018 www.quentinhuys.com qhuys@cantab.net

np = size(x,1);
if nargout==2; dodiff=1; else; dodiff=0;end


bmb  = exp(x(1));
bmf  = exp(x(2));
al = 1./(1+exp(-x(3)));
rep = x(4)*eye(2);

Q1 = zeros(2,1);
Q2 = zeros(2,2);
dQ1da= zeros(2,1);
dQ2da= zeros(2,2);
dQeda= zeros(2,1);
dQ1dl = zeros(2,1);
dQedl = zeros(2,1);
drep = eye(2);

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

if options.generatesurrogatedata==1
	dodiff=0;
	rewprob = D.rewprob; 
	trans = D.trans;
end

% if second-level states are coded as '2' and '3' change to '1' and '2' 
if any(D.S(2,:)==3); D.S(2,:) = D.S(2,:)-1; end

bb=20;
n=zeros(2);
for t=1:length(D.A);

	s=D.S(1,t); sp=D.S(2,t);
	a=D.A(1,t); ap=D.A(2,t);
	r=D.R(1,t);

	if ~isnan(a) & ~isnan(ap); 

		if n(1,1)+n(2,2) > n(1,2)+n(2,1)
			Tr = .3+.4*eye(2);
		else
			Tr = .3+.4*(1-eye(2));
		end

		pqm = bb*Q2;
		pqm = pqm-ones(2,1)*max(pqm);
		pqm = exp(pqm);
		pqm = pqm*diag(1./sum(pqm));
		Qd = Tr*sum(Q2.*pqm)';

		Qeff= bmb*Qd; 
		if t>1 & exist('a1old'); Qeff= Qeff + rep(:,a1old); end

		lpa = Qeff;
		lpa = lpa - max(lpa);
		lpa = lpa - log(sum(exp(lpa)));
		pa = exp(lpa);
		if options.generatesurrogatedata==1
			[a,sp] = simulateTwostep(pa,s,trans(t));
		else
			l = l + lpa(a);
		end

		lpap = bmf*Q2(:,sp);
		lpap = lpap - max(lpap);
		lpap = lpap - log(sum(exp(lpap)));
		pap = exp(lpap);
		if options.generatesurrogatedata==1
			[ap,r] = simulateTwostep(pap,sp,trans(t),rewprob(:,:,t));
		else
			l = l + lpap(ap);
		end

		de2 = r - Q2(ap,sp);

		if dodiff

			dl(1) = dl(1) + bmb*(Qd(a)-pa'*Qd);
			dl(2) = dl(2) + bmf*(Q2(ap,sp)-pap'*Q2(:,sp));
			dl(3) = dl(3) + bmf*(dQ2da(ap,sp) - pap'*dQ2da(:,sp));

			% grad wrt al(1)
			dpqmda = bb*pqm.*(dQ2da - ones(2,1)*sum(pqm.*dQ2da)); 
			dQdda  = Tr*sum(dQ2da.*pqm + Q2.*dpqmda)'; 
			dQ2da(ap,sp) = dQ2da(ap,sp) + al*(1-al)*de2 + al*-dQ2da(ap,sp);

			dl(3) = dl(3) + bmb*(dQdda(a) - pa'*dQdda); 
			% grad wrt rep 
			if t>1 & exist('a1old');
				dl(4) = dl(4) + ((a==a1old) - pa'*drep(:,a1old));
			end

		end

		Q2(ap,sp) = Q2(ap,sp) + al*de2;

		n(sp,a) = n(sp,a)+1;
		a1old = a; 

		if options.generatesurrogatedata==1
			dsurr.A(1,t)=a; dsurr.A(2,t)=ap;
			dsurr.S(1,t)=s; dsurr.S(2,t)=sp;
			dsurr.R(1,t)=r;
		end
	end

end
l=-l;
dl=-dl;
