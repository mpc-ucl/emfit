function [l,dl,dsurr] = llbgeq0(x,D,mu,nui,doprior,options);
% 
% Fit model to data from probabilistic reinforcement task (Pizzagalli et al. 2005
% Biological Psychiatry). 
% 
% This is the basic (Rescorla-Wagner) model. 
% 
% Mapping anhedonia onto reinforcement learning. A behavioural meta-analysis.
% Huys QJM, Pizzagalli DA, Bogdan R and Dayan P
% Biology of Mood & Anxiety Disorders (2013) 3:12 
%
% Quentin Huys 2016
% www.quentinhuys.com 

np = length(x);
dodiff = nargout==2; 

b = exp(x(1));					% reward sensitivity
g = exp(x(2));             % instruction sensitivity
eps = 1/(1+exp(-x(3)));    % learning rate 
q0 = x(4);                 % initial bias 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a;
r = D.r;
s = D.s;
I0 = [1 1; 0 0];
I = D.I;

if options.generatesurrogatedata==1
	a = NaN*zeros(size(a));
	r = NaN*zeros(size(a));
	dodiff=0;
end

dqde=zeros(2);
q   = q0*I0; 
n=zeros(1,2);

for t=1:length(a);
	if ~isnan(D.a(t))

	qe = b*q(:,s(t)) + g*I(:,s(t));
	q0 = max(qe);
	lpa = qe-q0 - log(sum(exp(qe-q0)));
	pa = exp(lpa);

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pa,s(t),D.prc,D.I);
	else
		l = l + lpa(a(t));
	end

	if dodiff
		dl(1) = dl(1) + b*(q(a(t),s(t)) - pa'*q(:,s(t)));
		dl(2) = dl(2) + g*(I(a(t),s(t)) - pa'*I(:,s(t)));
		dl(3) = dl(3) + b*(dqde(a(t),s(t)) - pa'*dqde(:,s(t)));
		dl(4) = dl(4) + b*(1-eps)^n(s(t))* (I0(a(t),s(t)) - pa'*I0(:,s(t))); 

		if I0(a(t),s(t)); n(s(t))=n(s(t))+1; end

		dqde(a(t),s(t)) = (1-eps)*(dqde(a(t),s(t)) + eps*( r(t) - q(a(t),s(t)) ));
	end
	q(a(t),s(t))    = q(a(t),s(t)) + eps*( r(t)-q(a(t),s(t)) );

	end
end

l = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

