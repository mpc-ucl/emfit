function [a,r] = genageqa(x,s,Z);
% 
% Generate data from probabilistic reinforcement task (Pizzagalli et al. 2005
% Biological Psychiatry) for model fitted by likelihood llageqa.m 
% 
% Mapping anhedonia onto reinforcement learning. A behavioural meta-analysis.
% Huys QJM, Pizzagalli DA, Bogdan R and Dayan P
% Biology of Mood & Anxiety Disorders (2013) 3:12 
%
% Quentin Huys 2016
% www.quentinhuys.com 

np = length(x);

c = exp(x(1));					% reward sensitivity
g = exp(x(2));             % instruction sensitivity
epc = 1/(1+exp(-x(3)));    % learning rate 
q0 = x(4);                 % initial bias 

I0 = sum(Z.I0,2)/2;
qa  = q0*I0;
n=0;
r = zeros(length(s),1);
a = zeros(length(s),1);

for t=1:length(a);

	if ~isnan(a(t)) & ~isnan(s(t))

		qe = g*Z.I(:,s(t)) + c*qa;
		pa = exp(qe);
		sp = sum(pa);
		if ~isinf(sp);
			pa = pa/sp;
		else
			[foo,i]=max(qe);pa(i)=1;pa(3-i)=0;
		end

		a(t) = sum(rand>cumsum([0 pa']));

		if Z.I(a(t),s(t))
			r(t) = rand<Z.prc(a(t));
		else
			r(t) = 0;
		end


		qa(a(t))    = qa(a(t))             + epc*( r(t) - qa(a(t)) );

	end

end
