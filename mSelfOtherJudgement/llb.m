function [l,dl,dsurr] = llb(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple bias model
% 
% Use this within emfit.m to tit RL type models to a group of subjects using
% EM. 
% 
% Quentin Huys 2021

dodiff=nargout==2;
np = size(x,1);

selfposbias = x(1:2);

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
wordval = D.wordval; 
avatval = D.avatval;

Q = zeros(2,1);

if options.generatesurrogatedata==1
	a = zeros(size(a));
	dodiff=0;
end


for t=1:length(a);
	if ~isnan(a(t));

		bl = 1+(t>48);

		wv = (-wordval(t)+3)/2;
		
		Q(1) = selfposbias(bl)*wordval(t); 
		Q(2) = 0;  
		q0 = Q; 

		l0 = q0-max(q0);
		l0 = l0 - log(sum(exp(l0)));
		p = exp(l0); 

		if options.generatesurrogatedata==1
			[a(t),r(t)] = generatera(p,wordval(t),avatval(t));
		end
		l = l+l0(a(t));

		if dodiff
			tmp = [wordval(t);0];  
			dl(bl)   = dl(bl) + tmp(a(t)) - p'*tmp;
		end

	else
		if options.generatesurrogatedata==1
			a(t) = NaN;
			r(t) = NaN;
		end
	end
end

l = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

