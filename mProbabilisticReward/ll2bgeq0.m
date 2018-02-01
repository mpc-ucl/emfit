function [l,dl] = ll2bgeq0(x,D,mu,nui,doprior,options);
% 
% Fit model to data from probabilistic reinforcement task (Pizzagalli et al. 2005
% Biological Psychiatry). 
% 
% The basic (Rescorla-Wagner) model is llbgeq0. This model additionally allows
% non-rewarded trials to be treated as punishing. 
% 
% Mapping anhedonia onto reinforcement learning. A behavioural meta-analysis.
% Huys QJM, Pizzagalli DA, Bogdan R and Dayan P
% Biology of Mood & Anxiety Disorders (2013) 3:12 
%
% Quentin Huys 2016
% www.quentinhuys.com 

np = length(x);
dodiff = nargout==2;

b = exp(x(1:2));				% reward and non-reward sensitivity. The latter is negative 
g = exp(x(3));					% instruction sensitivity
eps = 1/(1+exp(-x(4)));		% learning rate 
q0 = x(5);						% initial bias 


% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a;
r = D.r;
s = D.s;
I = D.I;

if options.generatesurrogatedata==1
	a = zeros(size(a));
	r = zeros(size(a));
	dodiff=0;
end

pa = zeros(2,1);
q   = q0*[1 1; 0 0];
dqde=zeros(2);
dqdb=zeros(2,2,2);
dqdq=[1 1;0 0];

for t=1:length(a);
	if ~isnan(a(t))

		er = b(2-r(t)) * (2*r(t)-1);

		qe = q(:,s(t)) + g*I(:,s(t));
		q0 = max(qe);
		lpa = qe-q0 - log(sum(exp(qe-q0)));
		pa = exp(lpa);

		if options.generatesurrogatedata==1
			[a(t),r(t)] = generatera(pa',s(t),Z);
		else
			l = l + lpa(a(t));
		end

		if dodiff

			for k=1:2
				dl(k) = dl(k) + dqdb(a(t),s(t),k) - pa'*dqdb(:,s(t),k);
				tmp = (r(t)==1 & k==1) | (r(t)==0 & k==2);
				dqdb(a(t),s(t),k) = (1-eps)*dqdb(a(t),s(t),k) + eps*er*tmp;
			end

			dl(3) = dl(3) + g*(I(a(t),s(t)) - pa'*I(:,s(t)));
			dl(4) = dl(4) + (dqde(a(t),s(t)) - pa'*dqde(:,s(t)));
			dl(5) = dl(5) + (dqdq(a(t),s(t)) - pa'*dqdq(:,s(t)));

			dqde(a(t),s(t)) = (1-eps)*(dqde(a(t),s(t)) + eps*( er - q(a(t),s(t)) ));
			dqdq(a(t),s(t)) = (1-eps)*dqdq(a(t),s(t));
		
		end

		q(a(t),s(t))    = q(a(t),s(t)) + eps*( er - q(a(t),s(t)) );

	end
end

l = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end


