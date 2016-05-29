function [a,r] = genbgeq0(x,s,Z);
% 
% Generate data from probabilistic reinforcement task (Pizzagalli et al. 2005
% Biological Psychiatry) for model fitted by likelihood llbgeq0.m 
% 
% Mapping anhedonia onto reinforcement learning. A behavioural meta-analysis.
% Huys QJM, Pizzagalli DA, Bogdan R and Dayan P
% Biology of Mood & Anxiety Disorders (2013) 3:12 
%
% Quentin Huys 2016
% www.quentinhuys.com 

np = length(x);

b = exp(x(1:2));				% reward and non-reward sensitivity. The latter is negative 
g = exp(x(3));					% instruction sensitivity
eps = 1/(1+exp(-x(4)));		% learning rate 
q0 = x(5);						% initial bias 

q   = q0*Z.I0; 
r = zeros(length(s),1);
a = zeros(length(s),1);

for t=1:length(a);

	er = b(2-r(t)) * (2*r(t)-1);

	qe = q(:,s(t)) + g*Z.I(:,s(t));
	pa = exp(qe);
	sp = sum(pa);
	if ~isinf(sp);
		pa = pa/sp;
	else
		[foo,i]=max(qe);pa(i)=1;pa(3-i)=0;
	end

	a(t) = sum(rand>cumsum([0 pa']));

	if Z.I(s(t),a(t))
		r(t) = rand<Z.prc(a(t));
	else
		r(t) = 0;
	end

	q(a(t),s(t))    = q(a(t),s(t)) + eps*( er - q(a(t),s(t)) );

end

