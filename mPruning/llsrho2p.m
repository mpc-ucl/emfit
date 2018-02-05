function [l,dl,a,r,Pa,Pa2] = llsrho2p(x,D,mu,nui,doprior,op);
% 
% [l,dl,a,r,Pa,Pa2] = llsrho2p(x,D,mu,nui,doprior,op);
% 
% log likelihood (l) and gradient (dl) of pruning model. Pruning model but 
% fitting separate reward sensitivities to separate outcomes. 
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
rho = x(1:4);							% sensitivities to outcomes 
pr = 1./(1+exp(-x(5:6)));			% prunign parameters. first is general (falt), second specific after large losses.

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[lp,dlp] = logGaussianPrior(x,mu,nui,doprior);

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

	clear dQdr;
	for kk=1:4
		dQdr(kk,:) = sum((ri==kk) .* prg .* prs);
	end
	dQ0dprg = sum(r    .* prs .* pr(1).^(pg-1).*pg * pr(1)*(1-pr(1)));
	dQ0dprs = sum(r    .* prg .* pr(2).^(ps-1).*ps * pr(2)*(1-pr(2)));

	tmp = Q0;
	l0 = max(tmp);
	la = tmp-l0 - log(sum(exp(tmp-l0)));
	l(k) = la(an);
	pa = exp(la);

	for kk=1:4
		dl(kk,k) = dQdr(kk,an) - dQdr(kk,:)*pa';
	end
	dl(5,k) = dQ0dprg(an) - dQ0dprg*pa';
	dl(6,k) = dQ0dprs(an) - dQ0dprs*pa';

	if nargout >=3; 
		Qp(k-Z.include+1,1:2^D.dn(k)) = Q0; 
		Pa(k) = pa(an);
		Pa2(k,i) = pa;
		dQdgg(k) = dQ0dprg(an);
		dQdgs(k) = dQ0dprs(an);
		Qa(k) = Q0(an);

		dlPdgg(k) = dQ0dprg(an) - dQ0dprg*pa';
		dlPdgs(k) = dQ0dprs(an) - dQ0dprs*pa';
		dPdgg(k) = pa(an)*(dQ0dprg(an) - dQ0dprg*pa');
		dPdgs(k) = pa(an)*(dQ0dprs(an) - dQ0dprs*pa');
	end
 
end
l  = -sum(l) - sum(lp);
dl = -sum(dl,2) - dlp;


