function [l,dl, dsurr] = llreweffchRel(x,D,mu,nui,doprior,options)

rch = D.rch;
ru  = D.ru;
ech = D.ech;
eu  = D.eu;
gr = D.gr;

np = length(x);
dodiff= nargout==2;


br1 = exp(x(1));
be1 = exp(x(2));
br2 = exp(x(3));
be2 = exp(x(4));

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

% l  = -1/2*(x'*x)/20;
% dl = -x/20; 
% l = 0; 
% dl = zeros(2,1);

if options.generatesurrogatedata==1
	a = NaN*zeros(length(rch));
	dodiff=0;
end

for t=1:length(rch)
    
    if gr(t)==1
    
        q = br1*[rch(t);ru(t)] - be1*[ech(t);eu(t)];
        
    elseif gr(t)==2
       
        q = br2*[rch(t);ru(t)] - be2*[ech(t);eu(t)];

    end
        
        q0 = max(q);
        l0 = q-q0 - log(sum(exp(q-q0)));
        p = exp(l0)';

    if options.generatesurrogatedata==1
		a(t) = rand<=p(1);   
    else
        l = l+l0(1);
    end
	

	if dodiff
        if gr(t)==1
            
            dqdbr1 =  br1*[rch(t);ru(t)]; 
            dqdbe1 = -be1*[ech(t);eu(t)]; 
            dl(1) = dl(1) + dqdbr1(1) - p*dqdbr1;
            dl(2) = dl(2) + dqdbe1(1) - p*dqdbe1;
            
        elseif gr(t)==2

            dqdbr2 =  br2*[rch(t);ru(t)]; 
            dqdbe2 = -be2*[ech(t);eu(t)]; 
            dl(3) = dl(3) + dqdbr2(1) - p*dqdbr2;
		    dl(4) = dl(4) + dqdbe2(1) - p*dqdbe2;
        
        end       
        
	end

end

l = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
end

end