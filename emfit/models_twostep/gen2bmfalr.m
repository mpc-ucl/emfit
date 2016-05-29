function [a,r,s,trans] = gen2bmfalr(x,Tr,trans,rewprob);
%
% Fit joint tree search and SARSA(lambda) model with separate learning rates and
% betas to two-step task. Note here the model-free weights at level one and two
% are allowed to differ, but there is only one learning rate. 
% 
% Quentin Huys, 2016 
% www.quentinhuys.com/code.html 
% www.quentinhuys.com/pub.html
% qhuys@cantab.net

np = size(x,1);
if nargout==2; dodiff=1; else; dodiff=0;end
T = length(rewprob); 

bmf  = exp(x(1));
b2   = exp(x(2));
al = 1./(1+exp(-x(3)));
la = 1./(1+exp(-x(4)));
rep = x(5)*eye(2);

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

	lpap = b2*Q2(:,s(2,t));
	lpap = lpap - max(lpap);
	lpap = lpap - log(sum(exp(lpap)));
	pap = exp(lpap);

	a(2,t) = sum(rand>cumsum([0 pap']));
	r(t) = rewprob(s(2,t),a(2,t),t)>rand;

	de1 = Q2(a(2,t),s(2,t))-Q1(a(1,t));
	de2 = r(t) - Q2(a(2,t),s(2,t));

	Q1(a(1,t)) = Q1(a(1,t)) + al*(de1 + la*de2);
	Q2(a(2,t),s(2,t))  = Q2(a(2,t),s(2,t))  + al*de2;

end
