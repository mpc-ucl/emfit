function [l,dl,dsurr] = llba2epb(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llba2epxb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with irreducible
% noise, constant bias, joint reward/loss sensitivity, and separate positively 
% constrained Pavlovian bias parameters for rewards and losses. 
% 
% Use this within emfit.m to tit RL type models to a group of subjects using EM. 
% 
% Guitart-Masip M, Huys QJM, Fuentemilla L, Dayan P, Duezel E and Dolan RJ
% (2012): Go and nogo learning in reward and punishment: Interactions between
% affect and effect.  Neuroimage 62(1):154-66 
%
% Quentin Huys 2011-2012 qhuys@gatsby.ucl.ac.uk

dodiff=nargout==2;
np = length(x);
beta 		= exp(x(1));				% sensitivity to rewards and losses
alfa 		= 1./(1+exp(-x(2)));		% learning rate
epsilon 	= exp(x(3:4));				% 'pavlovian' parameter. Weigth of Vcue into Qgo
bias 		= x(5);						% constant bias 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

V=zeros(1,4); 
Q=zeros(2,4); 

dQdb = zeros(2,4,2);
dQde = zeros(2,4);
dVdb = zeros(4,4);
dVde = zeros(4,1);

if options.generatesurrogatedata==1
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
	rho = sum(s(t)==[1 3]);

	q = Q(:,s(t)); 
	q(1) = q(1) + epsilon(2-rho) * V(s(t)) + bias;    % add Pavlovian effect 

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p = exp(la); 

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera_exact(p',s(t),sum(s(t)==s),options.session);
	end
	l = l + la(a(t));

	er = beta * r(t);

	if dodiff
		tmp = (dQdb(:,s(t)) + [epsilon(2-rho)*dVdb(s(t));0]);
		dl(1) = dl(1) +  (tmp(a(t)) - p'*tmp); 
		dQdb(a(t),s(t)) = (1-alfa)*dQdb(a(t),s(t)) + alfa*er;
		dVdb(     s(t)) = (1-alfa)*dVdb(     s(t)) + alfa*er;

		tmp = (dQde(:,s(t)) + [epsilon(2-rho)*dVde(s(t));0]);
		dl(2) = dl(2) + (tmp(a(t)) - p'*tmp); 
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);
		dVde(     s(t)) = (1-alfa)*dVde(     s(t)) + (er-V(     s(t)))*alfa*(1-alfa);

		dl(4-rho) = dl(4-rho) + epsilon(2-rho)*V(s(t))* ((a(t)==1)  - p(1));

		tmp = [1;0];
		dl(5) = dl(5) +  (tmp(a(t)) - p'*tmp); 
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t)));  
	V(s(t))      = V(s(t))      + alfa * (er - V(s(t)     ));

end
l  = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end


