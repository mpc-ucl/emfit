function [l,dl,dsurr] = ll2baqx(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = ll2baqx(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with irreducible
% noise, initial bias and separate reward and loss sensitivities
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
beta 		= exp(x(1:2));				% sensitivity to rewards and losses
alfa 		= 1./(1+exp(-x(3)));		% learning rate
qa			= x(4);						% initial value for Q go
g    		= 1/(1+exp(-x(5)));		% irreducible noise 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

V=zeros(1,4); 
Q=zeros(2,4); 
Q(1,:)=qa;

dQdb = zeros(2,4,2);
dQde = zeros(2,4);
dQdq = zeros(2,4); dQdq(1,:)=1;

if options.generatesurrogatedata==1
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
	rho = sum(s(t)==[1 3]);

	q = Q(:,s(t)); 

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t));
	end
	l = l + log(pg(a(t)));
	er = beta(2-rho) * r(t);

	if dodiff
		for k=1:2
			dl(k) = dl(k) + g*(p0(a(t)) * (dQdb(a(t),s(t),k) - p0'*dQdb(:,s(t),k))) / pg(a(t));
			tmp = (rho==1 & k==1) | (rho==0 & k==2);
			dQdb(a(t),s(t),k) = (1-alfa)*dQdb(a(t),s(t),k) + alfa*er*tmp;
		end

		dl(3) = dl(3) + g*(p0(a(t)) * (dQde(a(t),s(t)) - p0'*dQde(:,s(t)))) / pg(a(t));
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);

		dl(4) = dl(4) + g*(p0(a(t)) * (dQdq(a(t),s(t)) - p0'*dQdq(:,s(t)))) / pg(a(t));
		dQdq(a(t),s(t)) = (1-alfa)*dQdq(a(t),s(t)) ;

		dl(5) = dl(5) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t)));  

end
l  = -l;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

