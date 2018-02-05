function [l,dl,a,r,Pa,Pa2] = llsb(x,D,mu,nui,doprior,op);
% 
% [l,dl,a,r,Pa,Pa2] = llsb(x,D,mu,nui,doprior,op);
% 
% log likelihood (l) and gradient (dl) of full lookahead model 
% 
% Use this within emfit.m to tit RL type models to a group of subjects using EM. 
% 
% Huys QJM*, Eshel N*, O'Lions E, Sheridan L, Dayan P and Roiser JP (2012):
% Bonsai trees in your head: How the Pavlovian system sculpts goal-directed
% choices by pruning decision trees PLoS Comp Biol 8(3): e1002410 
%
% Huys QJM, Lally N, Faulkner P, Eshel N, Seifritz E, Gershman SJ, Dayan P and
% Roiser JP (2015): The interplay of approximate planning strategies PNAS,
% 112(10):3098-3103
%
% Quentin Huys 2011-2017 qhuys@cantab.net


np = length(x);
beta 		= exp(x(1));				% sensitivity to lookahead computations

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[lp,dlp] = logGaussianPrior(x,mu,nui,doprior);

Z = D.Z;

Q0 = squeeze(sum(Z.R));
for ss=1:6
	for dd=3:Z.dmax
		i=1:2^dd;
		tmp  = beta*Q0(i,dd-2,ss);
		l0 = max(tmp);
		la(i,dd-2,ss) = tmp-l0 - log(sum(exp(tmp-l0)));
	end
end

pa = exp(la);
an=D.an; % binary code for action sequence for all trials in one subject 
dn=D.dn; 
sn=D.sn; 

for k=1:length(an); 
	a=an(k);
	d=dn(k);
	s=sn(k);
	i=1:2^dn(k);

	l(k) = la(a,d-2,s);
	dl(k) = beta*Q0(a,d-2,s) - beta*(Q0(i,d-2,s)'*pa(i,d-2,s));

	Pa(k) = pa(a,d-2,s);
	Pa2(k,i) = pa(i,d-2,s);
end
l  = -sum(l) - sum(lp);
dl = -sum(dl,2) - dlp;

