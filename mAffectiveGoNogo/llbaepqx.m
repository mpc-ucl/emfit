function [l,dl,dsurr] = llbaepqx(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llbaepxb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with constant bias
% towards one action, irreducible noise and positive Pavlovian bias parameter
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
beta 		= exp(x(1));				% sensitivity to reward          
alfa 		= 1./(1+exp(-x(2)));		% learning rate
epsilon 	= exp(x(3));				% 'pavlovian' parameter. Weigth of Vcue into Qgo
q0			= x(4);						% initial bias 
g       	= 1/(1+exp(-x(5)));		% irreducible noise

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

V=zeros(1,4); 
Q=zeros(2,4); 
Q(1,:) = q0;

dQdb = zeros(2,4);
dQde = zeros(2,4);
dVdb = zeros(4,4);
dVde = zeros(4,1);
dQdq = zeros(2,4); dQdq(1,:)=1;

if options.generatesurrogatedata==1
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)

	q = Q(:,s(t)); 
	q(1) = q(1) + epsilon * V(s(t)) ;    % add Pavlovian effect 

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t));
	end
	l = l + log(pg(a(t)));

	er = beta * r(t);

	if dodiff
		tmp = (dQdb(:,s(t)) + [epsilon*dVdb(s(t));0]);
		dl(1) = dl(1) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQdb(a(t),s(t)) = (1-alfa)*dQdb(a(t),s(t)) + alfa*er;
		dVdb(     s(t)) = (1-alfa)*dVdb(     s(t)) + alfa*er;

		tmp = (dQde(:,s(t)) + [epsilon*dVde(s(t));0]);
		dl(2) = dl(2) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);
		dVde(     s(t)) = (1-alfa)*dVde(     s(t)) + (er-V(     s(t)))*alfa*(1-alfa);

		dl(3) = dl(3) + g*(p0(a(t))*epsilon*V(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));

		dl(4) = dl(4) + g*(p0(a(t)) * (dQdq(a(t),s(t)) - p0'*dQdq(:,s(t)))) / pg(a(t));
		dQdq(a(t),s(t)) = (1-alfa)*dQdq(a(t),s(t)) ;

		dl(5) = dl(5) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

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



