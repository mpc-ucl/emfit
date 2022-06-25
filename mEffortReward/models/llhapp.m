function [l,dl,dsurr] = llhapp(x,D,mu,nui,doprior,options)

hap=d.hap;
rch=D.rch;
ech=D.ech;
ru =D.ru ;
eu =D.eu ;

np = length(x);
dodiff= nargout==2;


hapmean = exp(x(1));
hapresp = exp(x(2));

%l = 0; 
%dl = zeros(size(x,1),1); 
l  = 1/2*(x'*x)/10; 
dl = x/10; 



for t=1:length(rch)
        
	if ~isnan(hap(t));%~=-1 

		happred(t) = hapmean + bt*t/80 + bavrr*avrate; 
		dhap = hap(t)-happred(t);
		l = l + dhap^2; 

		dhappreddhapmean = 1; 
		dhappreddbavrr = avrate; 
		dhappreddbt = t/80; 
		dhappreddalpha =  + bavrr*davratedalpha; 

		dl(1) = dl(1) - 2*dhap*dhappreddhapmean; 
		dl(2) = dl(2) - 2*dhap*dhappreddbavrr; 
		dl(3) = dl(3) - 2*dhap*dhappreddbt;
		dl(4) = dl(4) - 2*dhap*dhappreddalpha; 

	end

end

