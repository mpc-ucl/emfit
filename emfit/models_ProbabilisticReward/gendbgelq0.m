function [a,r] = genbgelq0(x,s,Z);
% 
% Generate data from probabilistic reinforcement task (Pizzagalli et al. 2005
% Biological Psychiatry) for model fitted by likelihood lldbgeq0.m 
% 
% Mapping anhedonia onto reinforcement learning. A behavioural meta-analysis.
% Huys QJM, Pizzagalli DA, Bogdan R and Dayan P
% Biology of Mood & Anxiety Disorders (2013) 3:12 
%
% Quentin Huys 2016
% www.quentinhuys.com 

np = length(x);

b = exp(x(1));					% reward sensitivity 
g = exp(x(2));             % instruction sensitivity
eps = 1/(1+exp(-x(3)));    % learning rate 
bl = 1/(1+exp(-x(4)));     % belief - mixes stimuli 
q0 = x(5);                 % initial bias 

q   = q0*Z.I0; 
r = zeros(length(s),1);
a = zeros(length(s),1);

for t=1:length(a);

	qe = b*(bl*q(:,s(t))+(1-bl)*q(:,3-s(t))) + g*Z.I(:,s(t));
	q0 = max(qe);
	lpa = qe-q0 - log(sum(exp(qe-q0)));
	pa = exp(lpa);

	a(t) = sum(rand>cumsum([0 pa']));

	if Z.I(s(t),a(t))
		r(t) = rand<Z.prc(a(t));
	else
		r(t) = 0;
	end

	q(a(t),s(t))    = q(a(t),s(t))   + bl    *eps*( r(t)-q(a(t),s(t)) );
	q(a(t),3-s(t))  = q(a(t),3-s(t)) + (1-bl)*eps*( r(t)-q(a(t),3-s(t)) );

end

