function [a,r,s,trans] = genbmfalr(x,Tr,trans,rewprob);
%
% Generate data from the model-free component fitted with llbmfalr.m of the
% model in Daw et al. 2011. Note this version assumes that betas at first and
% second levels take on the same value

np = size(x,1);
if nargout==2; dodiff=1; else; dodiff=0;end
T = length(rewprob); 

bmf  = exp(x(1));
al = 1./(1+exp(-x(2)));
la = 1./(1+exp(-x(3)));
rep = x(4)*eye(2);

Q1 = zeros(2,1);
Q2 = zeros(2,2);

bb=20;
s(1,1:T) = 1;
for t=1:T
	
	Qeff= bmf*Q1;
	if t>1; Qeff = Qeff + rep(:,a(1,t-1)); end

	lpa = Qeff;
	lpa = lpa - max(lpa);
	lpa = lpa - log(sum(exp(lpa)));
	pa = exp(lpa);

	a(1,t) = sum(rand>cumsum([0 pa']));
	if trans(t); if a(1,t)==1; s(2,t) = 1; else s(2,t) = 2; end
	else       ; if a(1,t)==1; s(2,t) = 2; else s(2,t) = 1; end
	end

	lpap = bmf*Q2(:,s(2,t));
	lpap = lpap - max(lpap);
	lpap = lpap - log(sum(exp(lpap)));
	pap = exp(lpap);

	a(2,t) = sum(rand>cumsum([0 pap']));
	r(t) = rewprob(s(2,t),a(2,t),t)>rand;

	de1 = Q2(a(2,t),s(2,t))-Q1(a(1,t));
	de2 = r(t) - Q2(a(2,t),s(2,t));

	Q1(a(1,t))         = Q1(a(1,t))         + al*(de1 + la*de2);
	Q2(a(2,t),s(2,t))  = Q2(a(2,t),s(2,t))  + al*de2;

end
