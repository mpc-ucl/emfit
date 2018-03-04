function [l,dl,dsurr] = llsp(x,D,mu,nui,doprior,options);
% 
% [l,dl,dsurr] = llsp(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of lookahead model with flat pruning. 
% 
% 
% Use this within emfit.m to fit RL type models to a group of subjects using EM. 
% 
% Huys QJM*, Eshel N*, O'Lions E, Sheridan L, Dayan P and Roiser JP (2012):
% Bonsai trees in your head: How the Pavlovian system sculpts goal-directed
% choices by pruning decision trees PLoS Comp Biol 8(3): e1002410 
%
% The neural basis of aversive Pavlovian guidance during planning Lally*
% N, Huys* QJM, Eshel N, Faulkner P, Dayan P and Roiser JP J. Neurosci. (2017)
% 37(42):10215-10229 https://quentinhuys.com/pub/LallyEa17-PruningfMRI.pdf
%
% Quentin Huys 2011-2017 qhuys@cantab.net


dodiff=nargout==2;
np = length(x);

beta 		= exp(x(1));				% sensitivity to lookahead computations
pr = 1./(1+exp(-x(2)));				% flat pruning probability 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

Z = D.Z;

% precompute action probabilities for all state-depth combinations 
for ss=1:6
	for dd=3:Z.dmax  
		i=1:2^dd;

		r = Z.R(1:dd,i,dd-2,ss);
		g = Z.D(1:dd,i,dd-2,ss);
      % sums over reward-matrix times gamma^depth
		Q0    (i,dd-2,ss) = sum(r .* pr.^ g); 
		dQ0dpr(i,dd-2,ss) = sum(r .* pr.^(g-1).*g * pr*(1-pr));

		tmp = beta*Q0(i,dd-2,ss);
		l0 = max(tmp);
		la(i,dd-2,ss) = tmp-l0 - log(sum(exp(tmp-l0)));
	end
end
pa = exp(la);

% load data 
an=D.an;  % binary code for action sequence for all trials in one subject 
dn=D.dn; 
sn=D.sn; 

if options.generatesurrogatedata==1
	an = zeros(size(an));
	dodiff=0;
end

for k=1:length(an);
	a=an(k); 
	d=dn(k);
	s=sn(k);
	i=1:2^dn(k);

	if options.generatesurrogatedata==1
		[ansurr(k)] = simulateActions(pa(i,d-2,s)');
		% Pa(k) = pa(a,d-2,s);
		% Pa2(k,i) = pa(i,d-2,s);
	else
		l = l + la(a,d-2,s);
	end

	if dodiff
		dl(1) = dl(1) + beta* Q0   (a,d-2,s) - beta*( Q0   (i,d-2,s)'*pa(i,d-2,s));
		dl(2) = dl(2) + beta*dQ0dpr(a,d-2,s) - beta*(dQ0dpr(i,d-2,s)'*pa(i,d-2,s));
	end
end
l  = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.an = an; 
end
