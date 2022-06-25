function [l,dl, dsurr] = llreweffchRelMB(x,D,mu,nui,doprior,options)

rch = D.rch;
ru  = D.ru;
ech = D.ech;
eu  = D.eu;
gr = D.gr;
vr = D.vr;
ve = D.ve;

np = length(x);
dodiff= nargout==2;


br1 = exp(x(1)+x(2));
be1 = exp(x(3)+x(4));
% br2 = exp(vr+x(1)-x(2));
% be2 = exp(ve+x(3)-x(4));

br2 = exp(x(1)-x(2));
be2 = exp(x(3)-x(4));

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
            
            dx1 =  br1*[rch(t);ru(t)]; 
            dx2 =  br1*[rch(t);ru(t)]; 
            dx3 = -be1*[ech(t);eu(t)]; 
            dx4 = -be1*[ech(t);eu(t)]; 
            
            dl(1) = dl(1) + dx1(1) - p*dx1;
            dl(2) = dl(2) + dx2(1) - p*dx2;
            dl(3) = dl(3) + dx3(1) - p*dx3;
            dl(4) = dl(4) + dx4(1) - p*dx4;
            
        elseif gr(t)==2

            dx1 =  br1*[rch(t);ru(t)]; 
            dx2 = -br1*[rch(t);ru(t)]; 
            dx3 = -be1*[ech(t);eu(t)]; 
            dx4 =  be1*[ech(t);eu(t)]; 
            
            dl(1) = dl(1) + dx1(1) - p*dx1;
            dl(2) = dl(2) + dx2(1) - p*dx2;
            dl(3) = dl(3) + dx3(1) - p*dx3;
            dl(4) = dl(4) + dx4(1) - p*dx4;
        
        end       
        
	end

end

l = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
end

end
