function [l,dl,dsurr] = llbmfalr(x,D,mu,nui,doprior,options);
% 
% [l,dl,dsurr] = llbmfalr(x,D,mu,nui,doprior,options);
%
% Fit simple SARSA(lambda) model with a single learning rate and
% beta to two-step task. 
% 
% X are the parameters for which to evaluate the likelihood D contains the data
% (see dataformat.txt or generate some data using generateExampleDataset.m).  MU
% is the prior mean and NUI is the prior inverse covariance matrix.  The DOPRIOR
% flag defines whether to apply a prior (1) or not (0). The variable
% OPTIONS.generatesurrogatedata defines whether to return the likelihood of data
% D (0) or whether to generate new surrogate data from the given parameters (1). 
% 
% Quentin Huys, 2018 www.quentinhuys.com qhuys@cantab.net

if nargout==2; dodiff=1; else; dodiff=0;end


bmf  = exp(x(1));
al = 1./(1+exp(-x(2)));
la = 1./(1+exp(-x(3)));
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
a1old=-1;
for t=1:length(D.A);

	
	s=D.S(1,t); sp=D.S(2,t);
	a=D.A(1,t); ap=D.A(2,t);
	r=D.R(1,t);

	if ~isnan(a) && ~isnan(ap); 

		if n(1,1)+n(2,2) > n(1,2)+n(2,1)
			Tr = .3+.4*eye(2);
		else
			Tr = .3+.4*(1-eye(2));
		end

		Qeff= bmf*Q1;
		if t>1 && a1old > 0;
			Qeff= Qeff + rep(:,a1old); 
		end

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


		de1 = Q2(ap,sp)-Q1(a);
		de2 = r - Q2(ap,sp);


		if dodiff

			dl(1) = dl(1) + bmf*(Q1(a)-pa'*Q1 + Q2(ap,sp)-pap'*Q2(:,sp));

			dl(2) = dl(2) + bmf*(dQ1da(a) - pa'*dQ1da + dQ2da(ap,sp) - pap'*dQ2da(:,sp));

			dl(3) = dl(3) + dQedl(a)  - pa'*dQedl;

			% grad wrt al(1)
			dQ1da(a) = dQ1da(a) + al*(1-al)*(de1+la*de2) + al*(dQ2da(ap,sp)-dQ1da(a) + la*-dQ2da(ap,sp)); 
			dQ2da(ap,sp) = dQ2da(ap,sp) + al*(1-al)*de2 + al*-dQ2da(ap,sp);

			% grad wrt la
			dQ1dl(a) = dQ1dl(a) + al*(-dQ1dl(a) + la*(1-la)*de2);
			dQedl = bmf*dQ1dl;

			% grad wrt rep 
			if t>1 && a1old > 0;
				dl(4) = dl(4) + ((a==a1old) - pa'*drep(:,a1old));
			end

		end

		Q1(a)     = Q1(a)     + al*(de1 + la*de2);
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
