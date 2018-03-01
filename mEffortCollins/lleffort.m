function [l,dl,dsurr] = llconstant(x,D,mu,nui,doprior,options);
% 
% [l,dl,dsurr] = llconstant(x,D,mu,nui,doprior,options);
% 
% Fit model to data from effort task by Gold, J. M.; Strauss, G. P.; Waltz, J.
% A.; Robinson, B. M.; Brown, J. K. & Frank, M. J. Negative symptoms of
% schizophrenia are associated with abnormal effort-cost computations. Biol
% Psychiatry, 2013, 74, 130-136
% 
% This is the most basic model which assumes participants entirely disregard the
% reward information and only focus on the effort. It contains one parameter
% theta measuring the effort sensitivity. 
% 
% Quentin Huys 2018 www.quentinhuys.com 

np = length(x);
dodiff= nargout==2;

betaeff = exp(x(1));	% beta for efforts 


% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.rew; 

effortCostLo = D.effortCostLo;%  set to - 0.2 in original task; 
effortCostHi = D.effortCostHi;%  set to - 1 in original task; 
effortCost = [effortCostLo effortCostHi]';

if options.generatesurrogatedata==1
	a = NaN*zeros(size(a));
	dodiff=0;
end

for t=1:D.Nch
	qe = betaeff*effortCost;

	% softmax avoiding numerical issues - subtract maximum
	q0 = max(qe);
	lpa = qe-q0 - log(sum(exp(qe-q0)));
	pa = exp(lpa);

	if options.generatesurrogatedata==1
		a(t) = simulateEffort(pa);
	else
		l = l + lpa(a(t));
	end

	if dodiff
		dl(1) = dl(1) + betaeff*(effortCost(a(t)) - pa'*effortCost);
	end

end
l = -l ; 
dl = - dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
end
