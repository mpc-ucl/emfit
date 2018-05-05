function [l,dl,dsurr] = llageqa(x,D,mu,nui,doprior,options);
% 
% Fit model to data from probabilistic reinforcement task (Pizzagalli et al. 2005
% Biological Psychiatry). 
% 
% The basic (Rescorla-Wagner) model is llbgeq0. This model disregards the
% different stimuli during learning, i.e. assumes that subjects learn only about
% how actions are rewarded, not how action-stimulus combinations are rewarded. 
% 
% Mapping anhedonia onto reinforcement learning. A behavioural meta-analysis.
% Huys QJM, Pizzagalli DA, Bogdan R and Dayan P
% Biology of Mood & Anxiety Disorders (2013) 3:12 
%
% Quentin Huys 2016
% www.quentinhuys.com 

np = length(x);
dodiff= nargout==2;

c = exp(x(1));					% reward sensitivity
g = exp(x(2));             % instruction sensitivity
epc = 1/(1+exp(-x(3)));    % learning rate 
q0 = x(4);                 % initial bias 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a;
r = D.r;
s = D.s;
I = D.I;
I0 = [1;0];

if options.generatesurrogatedata==1
	a = zeros(size(a));
	r = zeros(size(a));
	dodiff=0;
end

dqade=zeros(2,1);
qa  = q0*I0;
n=0;

for t=1:length(a);
	if ~isnan(a(t))

	qe = g*I(:,s(t)) + c*qa;
	q0 = max(qe);
	lpa = qe-q0 - log(sum(exp(qe-q0)));
	pa = exp(lpa);

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pa,s(t),D.prc,D.I);
	else
		l = l + lpa(a(t));
	end

	if dodiff
		dl(1) = dl(1) + c*(qa(a(t))        - pa'*qa);
		dl(2) = dl(2) + g*(I(a(t),s(t))  - pa'*I(:,s(t)));
		dl(3) = dl(3) + c*(dqade(a(t))  - pa'*dqade);

		dl(4) = dl(4) + c*(1-epc)^n* (I0(a(t)) - pa'*I0); 
		if I0(a(t)); n=n+1;end
		dqade(a(t)) = (1-epc)*(dqade(a(t)) + epc*( r(t) - qa(a(t))));
	end

	qa(a(t))    = qa(a(t))             + epc*( r(t) - qa(a(t)) );

	end
end

l = -l;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

