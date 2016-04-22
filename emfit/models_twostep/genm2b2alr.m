function [a,r,s,trans] = genm2b2alr(x,rewprob);
%
% Generate data from joint tree search and SARSA(lambda) model (fitted with
% llm2b2alr.m) with separate learning rates and betas to two-step task, as in Daw
% et al. 2011

np = size(x,1);
if nargout==2; dodiff=1; else; dodiff=0;end
T = length(rewprob); 

b  = exp(x(1:2));
al = 1./(1+exp(-x(3:4)));
la = 1./(1+exp(-x(5)));
w  = 1./(1+exp(-x(6)));
rep = x(7)*eye(2);

Q1 = zeros(2,1);
Q2 = zeros(2,2);

bb=20;
Tr = .3+.4*eye(2);
s(1,1:T) = 1;
n = zeros(2);
trans(1,1:T) = rand(1,T)>.3; 
for t=1:T
	
	if n(1,1)+n(2,2) > n(1,2)+n(2,1)
		Tr = .3+.4*eye(2);
	else
		Tr = .3+.4*(1-eye(2));
	end

	pqm = bb*Q2;
	pqm = pqm-ones(2,1)*max(pqm);
	pqm = exp(pqm);
	pqm = pqm*diag(1./sum(pqm));
	Qd = (sum(Q2.*pqm)*Tr)';
	%Qd = (max(Q2)*Tr)';

	if t>1
		Qeff= w*Qd + (1-w)*Q1 + rep(:,a(1,t-1));
	else
		Qeff= w*Qd + (1-w)*Q1;
	end
	lpa = b(1)*Qeff;
	lpa = lpa - max(lpa);
	lpa = lpa - log(sum(exp(lpa)));
	pa = exp(lpa);

	a(1,t) = sum(rand>cumsum([0 pa']));
	if trans(t); if a(1,t)==1; sp = 1; else sp = 2; end
	else       ; if a(1,t)==1; sp = 2; else sp = 1; end
	end

	lpap = b(2)*Q2(:,sp);
	lpap = lpap - max(lpap);
	lpap = lpap - log(sum(exp(lpap)));
	pap = exp(lpap);

	ap = sum(rand>cumsum([0 pap']));
	s(2,t) = sp;
	a(2,t) = ap;
	r(t) = rewprob(s(2,t),a(2,t),t)>rand;

	de1 = Q2(ap,sp)-Q1(a(1,t));
	de2 = r(t) - Q2(ap,sp);

	Q1(a(1,t)) = Q1(a(1,t)) + al(1)*(de1 + la*de2);
	Q2(ap,sp)  = Q2(ap,sp)  + al(2)*de2;

	n(sp,a(1,t)) = n(sp,a(1,t))+1;

end
