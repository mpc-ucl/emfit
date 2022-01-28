function [l,dl,dsurr] = llreweffchfat(x,D,mu,nui,doprior,options)

rch = D.rch;
ru  = D.ru;
ech = D.ech;
eu  = D.eu;

np = length(x);
dodiff= nargout==2;


br = exp(x(1));
be = exp(x(2));

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

% l  = -1/2*(x'*x)/20;
% dl = -x/20; 
%l = 0; 
%dl = zeros(2,1);

if options.generatesurrogatedata==1
	a = NaN*zeros(length(rch));
	dodiff=0;
end

for t=1:length(rch)

	q = br*[rch(t);ru(t)] - be*(t/80)*[ech(t);eu(t)];
	q0 = max(q);
	l0 = q-q0 - log(sum(exp(q-q0)));
	p = exp(l0)';
    
    if options.generatesurrogatedata==1
		a(t) = rand<=p(1);   
    else
        l = l+l0(1);
    end
	

	if dodiff
		dqdbr =  br*[rch(t);ru(t)]; 
        dqdbe = -be*(t/80)*[ech(t);eu(t)]; 
        
        dl(1) = dl(1) + dqdbr(1) - p*dqdbr;
		dl(2) = dl(2) + dqdbe(1) - p*dqdbe;
	end

end

l = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
end

end
