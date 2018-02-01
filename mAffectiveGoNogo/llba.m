function [l,dl,dsurr] = llba(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llba(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model 
% 
% Use this within emfit.m to tit RL type models to a group of subjects using
% EM. 
% 
% Guitart-Masip M, Huys QJM, Fuentemilla L, Dayan P, Duezel E and Dolan RJ
% (2012): Go and nogo learning in reward and punishment: Interactions between
% affect and effect.  Neuroimage 62(1):154-66 
%
% Quentin Huys 2018 qhuys@cantab.net

dodiff=nargout==2;
np = length(x);
beta 		= exp(x(1));				% sensitivity to reward          
alfa 		= 1./(1+exp(-x(2)));		% learning rate

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

V=zeros(1,4); 
Q=zeros(2,4); 

dQdb = zeros(2,4);
dQde = zeros(2,4);

if options.generatesurrogatedata==1
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)

	q = Q(:,s(t)); 

	l0 = max(q);
	la = q - l0 - log(sum(exp(q-l0)));
	pa = exp(la);

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pa',s(t));
	end
	l = l + la(a(t));

	er = beta * r(t);

	if dodiff
		dl(1) = dl(1) + dQdb(a(t),s(t)) - pa'*dQdb(:,s(t));
		dl(2) = dl(2) + dQde(a(t),s(t)) - pa'*dQde(:,s(t));

		dQdb(a(t),s(t)) = (1-alfa)*dQdb(a(t),s(t)) + alfa*er;
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t)));  

end
l  = -l ;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

