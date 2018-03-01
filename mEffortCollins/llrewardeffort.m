function [l,dl,dsurr] = llreweffscaling(x,D,mu,nui,doprior,options)
% 
% [l,dl,dsurr] = llreweffscaling(x,D,mu,nui,doprior,options)
% 
% Fit model to data from effort task by Anne Collins. 
% 
% This is the most basic model which assumes participants entirely disregard the
% reward information and only focus on the effort. It contains one parameter
% theta measuring the effort sensitivity. 
% 
% Quentin Huys 2018
% www.quentinhuys.com 

np = length(x);
dodiff= nargout==2;

betarew = exp(x(1));	% beta for rewards 
betaeff = exp(x(2));	% beta for efforts 

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
	qe = [betarew*[1;r(t)] + betaeff*effortCost];

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
		dqdbr = betarew*[1 r(t)]';
		dqdbe = betaeff*effortCost;
		dl(1) = dl(1) + dqdbr(a(t)) - pa'*dqdbr;
		dl(2) = dl(2) + dqdbe(a(t)) - pa'*dqdbe;
	end

end
l = -l; 
dl = - dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
end
