function [l,dl] = ll2bmf2alr(x,D,mu,nui,doprior);
%
% Fit joint tree search and SARSA(lambda) model with separate learning rates and
% betas to two-step task. Note here both the learning rates and the model-free
% weights at level one and two are allowed to differ. 
%  
% Quentin Huys, 2016 
% www.quentinhuys.com/code.html 
% www.quentinhuys.com/pub.html
% qhuys@cantab.net

np = size(x,1);
if nargout==2; dodiff=1; else; dodiff=0;end


bmf  = exp(x(1));
b2   = exp(x(2));
al = 1./(1+exp(-x(3:4)));
la = 1./(1+exp(-x(5)));
rep = x(6)*eye(2);

Q1 = zeros(2,1);
Q2 = zeros(2,2);
dQ1da1= zeros(2,1);
dQ1da2= zeros(2,1);
dQeda1= zeros(2,1);
dQeda2= zeros(2,1);
dQ2da2= zeros(2,2);
dQ1dl = zeros(2,1);
dQedl = zeros(2,1);
dQedw = zeros(2,1);
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

		if t>1 & exist('a1old');
			Qeff= bmf*Q1 + rep(:,a1old);
		else
			Qeff= bmf*Q1;
		end
		lpa = Qeff;
		lpa = lpa - max(lpa);
		lpa = lpa - log(sum(exp(lpa)));
		l = l + lpa(a);
		pa = exp(lpa);

		lpap = b2*Q2(:,sp);
		lpap = lpap - max(lpap);
		lpap = lpap - log(sum(exp(lpap)));
		l = l + lpap(ap);
		pap = exp(lpap);

		de1 = Q2(ap,sp)-Q1(a);
		de2 = r - Q2(ap,sp);

		if dodiff
			dl(1) = dl(1) + bmf*(Q1(a) - pa'*Q1);
			dl(2) = dl(2) + b2*(Q2(ap,sp) - pap'*Q2(:,sp));
			dl(3) = dl(3) + dQeda1(a) - pa'*dQeda1;
			dl(4) = dl(4) + b2*(dQ2da2(ap,sp)- pap'*dQ2da2(:,sp)) + bmf*(dQ1da2(a) - pa'*dQ1da2);
			dl(5) = dl(5) + (dQedl(a)  - pa'*dQedl);

			% grad wrt al(1)
			dQ1da1(a) = dQ1da1(a) + al(1)*(1-al(1))*(de1+la*de2) + al(1)*-dQ1da1(a);
			dQeda1 = bmf*dQ1da1;

			% grad wrt al(2)
			dQ1da2(a) = dQ1da2(a) + al(1)*(dQ2da2(ap,sp)-dQ1da2(a) + la*-dQ2da2(ap,sp));
			dQ2da2(ap,sp) = dQ2da2(ap,sp) + al(2)*(1-al(2))*de2 + al(2)*-dQ2da2(ap,sp);

			% grad wrt la
			dQ1dl(a) = dQ1dl(a) + al(1)*la*(1-la)*de2+ al(1)*-dQ1dl(a);
			dQedl = bmf*dQ1dl;

			% grad wrt rep 
			if t>1 & exist('a1old');
				dl(6) = dl(6) + ((a==a1old) - pa'*drep(:,a1old));
			end

		end

		Q1(a)     = Q1(a)     + al(1)*(de1 + la*de2);
		Q2(ap,sp) = Q2(ap,sp) + al(2)*de2;

		n(sp,a) = n(sp,a)+1;
		a1old = a; 
	end

end
l=-l;
dl=-dl;
