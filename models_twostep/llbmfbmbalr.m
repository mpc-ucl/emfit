function [l,dl] = llbmfbmbalr(x,D,mu,nui,doprior);
%
% Fit joint tree search and SARSA(lambda) model with separate learning rates and
% betas to two-step task. This version is reparametrised such that there are two
% separate weights for model-based and model-free components, rather than a
% weight explicitly trading off the two components as in Daw et al. 2011. Note
% here the model-free weights at level one and two are assumed to be equal, and
% there is only one learning rate. 
%
% Quentin Huys, 2015 
% www.quentinhuys.com/code.html 
% www.quentinhuys.com/pub.html
% qhuys@cantab.net

np = size(x,1);
if nargout==2; dodiff=1; else; dodiff=0;end


bmb  = exp(x(1));
bmf  = exp(x(2));
al = 1./(1+exp(-x(3)));
la = 1./(1+exp(-x(4)));
rep = x(5)*eye(2);

Q1 = zeros(2,1);
Q2 = zeros(2,2);
dQ1da= zeros(2,1);
dQ2da= zeros(2,2);
dQeda= zeros(2,1);
dQ1dl = zeros(2,1);
dQedl = zeros(2,1);
drep = eye(2);

if doprior;
	l = -1/2*(x-mu)'*nui*(x-mu) - np/2*log(2*pi) - 1/2*log(1/det(nui));
	dl = - nui*(x-mu);
else;
	l=0;
	dl=zeros(np,1);
end

bb=20;
n=zeros(2);
for t=1:length(D.A);

	
	s=D.S(1,t); sp=D.S(2,t)-1;
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

		Qeff= bmb*Qd + bmf*Q1;
		if t>1 & exist('a1old'); Qeff= Qeff + rep(:,a1old); end

		lpa = Qeff;
		lpa = lpa - max(lpa);
		lpa = lpa - log(sum(exp(lpa)));
		l = l + lpa(a);
		pa = exp(lpa);

		lpap = bmf*Q2(:,sp);
		lpap = lpap - max(lpap);
		lpap = lpap - log(sum(exp(lpap)));
		l = l + lpap(ap);
		pap = exp(lpap);

		de1 = Q2(ap,sp)-Q1(a);
		de2 = r - Q2(ap,sp);

		if dodiff
			dl(1) = dl(1) + bmb*(Qd(a)-pa'*Qd);
			dl(2) = dl(2) + bmf*(Q1(a)-pa'*Q1 + Q2(ap,sp)-pap'*Q2(:,sp));

			dl(3) = dl(3) + bmf*(dQ1da(a) - pa'*dQ1da + dQ2da(ap,sp) - pap'*dQ2da(:,sp));

			dl(4) = dl(4) + dQedl(a)  - pa'*dQedl;

			% grad wrt al(1)
			dpqmda = bb*pqm.*(dQ2da - ones(2,1)*sum(pqm.*dQ2da)); 
			dQdda  = Tr*sum(dQ2da.*pqm + Q2.*dpqmda)'; 
			dQ1da(a) = dQ1da(a) + al*(1-al)*(de1+la*de2) + al*(dQ2da(ap,sp)-dQ1da(a) + la*-dQ2da(ap,sp)); 
			dQ2da(ap,sp) = dQ2da(ap,sp) + al*(1-al)*de2 + al*-dQ2da(ap,sp);

			% grad wrt la
			dQ1dl(a) = dQ1dl(a) + al*(-dQ1dl(a) + la*(1-la)*de2);
			dQedl = bmf*dQ1dl;

			dl(3) = dl(3) + bmb*(dQdda(a) - pa'*dQdda); 

			% grad wrt rep 
			if t>1 & exist('a1old');
				dl(5) = dl(5) + ((a==a1old) - pa'*drep(:,a1old));
			end

		end

		Q1(a)     = Q1(a)     + al*(de1 + la*de2);
		Q2(ap,sp) = Q2(ap,sp) + al*de2;

		n(sp,a) = n(sp,a)+1;
		a1old = a; 
	end

end
l=-l;
dl=-dl;
