function [l,dl,dsurr] = llbaqx(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llbaqx(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with irreducible noise
% and with initial bias. 
% 
% Use this within emfit.m to tit RL type models to a group of subjects using
% EM. 
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
qa   		= x(3);						% initial value for Q go
g    		= 1/(1+exp(-x(4)));		% irreducible noise 

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

	q = Q(:,s(t)); 

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
		dl(1) = dl(1) + g*(p0(a(t)) * (dQdb(a(t),s(t)) - p0'*dQdb(:,s(t)))) / pg(a(t));
		dl(2) = dl(2) + g*(p0(a(t)) * (dQde(a(t),s(t)) - p0'*dQde(:,s(t)))) / pg(a(t));
		dl(3) = dl(3) + g*(p0(a(t)) * (dQdq(a(t),s(t)) - p0'*dQdq(:,s(t)))) / pg(a(t));
		dl(4) = dl(4) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

		dQdb(a(t),s(t)) = (1-alfa)*dQdb(a(t),s(t)) + alfa*er;
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);
		dQdq(a(t),s(t)) = (1-alfa)*dQdq(a(t),s(t)) ;
	end


	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t)));  

end
l  = -l;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

