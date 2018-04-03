function [l,dl,dsurr] = llsrho2p(x,D,mu,nui,doprior,options);
% 
% [l,dl,dsurr] = llsrho2p(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of pruning model with separate reward
% sensitivities to separate outcomes. 
% 
% Use this within emfit.m to tit RL type models to a group of subjects using EM. 
% 
% Huys QJM*, Eshel N*, O'Lions E, Sheridan L, Dayan P and Roiser JP (2012):
% Bonsai trees in your head: How the Pavlovian system sculpts goal-directed
% choices by pruning decision trees PLoS Comp Biol 8(3): e1002410 
%
% The neural basis of aversive Pavlovian guidance during planning Lally* N,
% Huys* QJM, Eshel N, Faulkner P, Dayan P and Roiser JP J. Neurosci. (2017)
% 37(42):10215-10229 https://quentinhuys.com/pub/LallyEa17-PruningfMRI.pdf
%
% Quentin Huys 2018 qhuys@cantab.net

dodiff=nargout==2;
np = length(x);

rho = x(1:4);							% sensitivities to outcomes 
pr = 1./(1+exp(-x(5:6)));			% pruning parameters. first is general (flat), second specific after large losses.

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

Z = D.Z;

for k=1:length(D.an);

	dn=D.dn(k); 
	sn=D.sn(k); 
	an=D.an(k);
	i=1:2^dn;

	ri = Z.Ri(1:dn,i,dn-2,sn);
	pg = Z.Pg(1:dn,i,dn-2,sn);
	ps = Z.Ps(1:dn,i,dn-2,sn);

	r = rho(ri);
	prg = pr(1).^pg;  % x(5) is specific continuation probability 
	prs = pr(2).^ps;	% x(6) is specific continuation probability 

	Q0      = sum(r    .* prg .* prs ); 

	if dodiff
		clear dQdr;
		for kk=1:4
			dQdr(kk,:) = sum((ri==kk) .* prg .* prs);
		end
		dQ0dprg = sum(r    .* prs .* pr(1).^(pg-1).*pg * pr(1)*(1-pr(1)));
		dQ0dprs = sum(r    .* prg .* pr(2).^(ps-1).*ps * pr(2)*(1-pr(2)));
	end

	tmp = Q0;
	l0 = max(tmp);
	la = tmp-l0 - log(sum(exp(tmp-l0)));
	pa = exp(la);

	if options.generatesurrogatedata==1
		[ansurr(k,1)] = simulateActions(pa);
		% Pa(k) = pa(an);
		% Pa2(k,i) = pa;
	else
		l = l + la(an);
	end

	if dodiff
		for kk=1:4
			dl(kk) = dl(kk) + dQdr(kk,an) - dQdr(kk,:)*pa';
		end
		dl(5) = dl(5) + dQ0dprg(an) - dQ0dprg*pa';
		dl(6) = dl(6) + dQ0dprs(an) - dQ0dprs*pa';
	end
 
end
l  = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.an = ansurr; 
end
