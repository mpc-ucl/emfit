function [l,dl] = lldbgelq0(x,D,mu,nui,doprior,options);
% 
% Fit model to data from probabilistic reinforcement task (Pizzagalli et al. 2005
% Biological Psychiatry). 
% 
% The basic (Rescorla-Wagner) model is llbgeq0. This model 
%
%	- additionally 'mixes' the stimuli to account for the fact that participants
%	  are unsure about which stimulus was presented. 
%	- additionally does counterfactual updates, i.e. learns about both stimuli on
%	  each trial, weighted by how likely they were to have been presented. 
% 
% Mapping anhedonia onto reinforcement learning. A behavioural meta-analysis.
% Huys QJM, Pizzagalli DA, Bogdan R and Dayan P 
% Biology of Mood & Anxiety Disorders (2013) 3:12 
%
% Quentin Huys 2016 
% www.quentinhuys.com 


np = length(x);
dodiff=nargout==2;

b = exp(x(1));					% reward sensitivity 
g = exp(x(2));             % instruction sensitivity
eps = 1/(1+exp(-x(3)));    % learning rate 
bl = 1/(1+exp(-x(4)));     % belief - mixes stimuli 
q0 = x(5);                 % initial bias 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a;
r = D.r;
s = D.s;
I0 = [1 1; 0 0];
I = D.I;

if options.generatesurrogatedata==1
	a = zeros(size(a));
	r = zeros(size(a));
	dodiff=0;
end

dqde=zeros(2);
dqdbl=zeros(2);
dqdq=I0;
q   = q0*I0; 

for t=1:length(a);
	if ~isnan(a(t))

	qe = b*(bl*q(:,s(t))+(1-bl)*q(:,3-s(t))) + g*I(:,s(t));
	q0 = max(qe);
	lpa = qe-q0 - log(sum(exp(qe-q0)));
	pa = exp(lpa);

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pa',s(t),Z);
	else
		l = l + lpa(a(t));
	end

	if dodiff;
		dl(1) = dl(1) + b*(bl*q(a(t),s(t))+(1-bl)*q(a(t),3-s(t)) - pa'*(bl*q(:,s(t))+(1-bl)*q(:,3-s(t))));
		dl(2) = dl(2) + g*(I(a(t),s(t)) - pa'*I(:,s(t)));
		dl(3) = dl(3) + b*(bl*dqde(a(t),s(t)) + (1-bl)*dqde(a(t),3-s(t)) - pa'*(bl*dqde(:,s(t))+(1-bl)*dqde(:,3-s(t))));

		tmp = bl*(1-bl)*(q(:,s(t))-q(:,3-s(t))) + bl*dqdbl(:,s(t))+(1-bl)*dqdbl(:,3-s(t));
		dl(4) = dl(4) + b*(tmp(a(t)) - pa'*tmp);

		dl(5) = dl(5) + b*(bl*dqdq(a(t),s(t)) + (1-bl)*dqdq(a(t),3-s(t)) - pa'*(bl*dqdq(:,s(t))+(1-bl)*dqdq(:,3-s(t))));

		dqdbl(a(t),  s(t)) = dqdbl(a(t),  s(t)) + bl*(1-bl)*eps*(r(t)-q(a(t),  s(t))) +    bl *eps*-dqdbl(a(t),  s(t));
		dqdbl(a(t),3-s(t)) = dqdbl(a(t),3-s(t)) - bl*(1-bl)*eps*(r(t)-q(a(t),3-s(t))) + (1-bl)*eps*-dqdbl(a(t),3-s(t));

		dqde(a(t),  s(t)) = (1-   bl *eps)*dqde(a(t),s(t))      + bl *eps*(1-eps)*( r(t) - q(a(t),  s(t)) );
		dqde(a(t),3-s(t)) = (1-(1-bl)*eps)*dqde(a(t),3-s(t)) + (1-bl)*eps*(1-eps)*( r(t) - q(a(t),3-s(t)) );
		dqdq(a(t),  s(t)) = (1-   bl *eps)*dqdq(a(t),s(t)); 
		dqdq(a(t),3-s(t)) = (1-(1-bl)*eps)*dqdq(a(t),3-s(t)); 
	end
	
	q(a(t),s(t))    = q(a(t),s(t))   + bl    *eps*( r(t)-q(a(t),s(t)) );
	q(a(t),3-s(t))  = q(a(t),3-s(t)) + (1-bl)*eps*( r(t)-q(a(t),3-s(t)) );

	end

end

l = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

