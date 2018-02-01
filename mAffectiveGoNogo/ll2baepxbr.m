function [l,dl] = ll2baepxb(x, a, r, s, Z,doprior)
dodiff=nargout==2;
np = length(x);
beta 	  = exp(x(1:2));            % sensitivity to reward          
alfa 	  = 1./(1+exp(-x(3)));     % learning rate
epsilon = exp(x(4));              % 'pavlovian' parameter. Weigth of Vcue into Qgo
g       = 1/(1+exp(-x(5)));
bias	  = x(6);
rep     = x(7);

if doprior
	l = -1/2 * (x-Z.mu)'*Z.nui*(x-Z.mu) - np/2*log(2*pi) - 1/2*log(1/det(Z.nui)); %
	dl = - Z.nui*(x-Z.mu) ;
else
	l = 0;
	dl = zeros(np,1);
end

repm = eye(2);
pa=zeros(2,1);
V=zeros(1,4); 
Q=zeros(2,4); 

dQdb = zeros(2,4,2);
dQde = zeros(2,4);
dVdb = zeros(4,2);
dVde = zeros(4,1);

for t=1:length(a)
	rho = sum(s(t)==[1 3]);
	er = beta(2-rho) * r(t);

	q = Q(:,s(t)); 
	q(1) = q(1) + epsilon * V(s(t)) + bias;    % add Pavlovian effect 
	if t>1
		q = q + x(7)*repm(:,a(t-1));
	end

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	l = l + log(pg(a(t)));

	if dodiff
		for k=1:2
			tmp = (dQdb(:,s(t),k) + [epsilon*dVdb(s(t),k);0]);
			dl(k) = dl(k) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
			tmp = (rho==1 & k==1) | (rho==0 & k==2);
			dQdb(a(t),s(t),k) = (1-alfa)*dQdb(a(t),s(t),k) + alfa*er * tmp;
			dVdb(     s(t),k) = (1-alfa)*dVdb(     s(t),k) + alfa*er * tmp;
		end

		tmp = (dQde(:,s(t)) + [epsilon*dVde(s(t));0]);
		dl(3) = dl(3) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);
		dVde(     s(t)) = (1-alfa)*dVde(     s(t)) + (er-V(     s(t)))*alfa*(1-alfa);

		dl(4) = dl(4) + g*(p0(a(t))*epsilon*V(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));

		dl(5) = dl(5) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(6) = dl(6) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));

		if t>1
			tmp = repm(:,a(t-1));
			dl(7) = dl(7) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		end
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t)));  
	V(s(t))      = V(s(t))      + alfa * (er - V(s(t)     ));

end
l  = -l; 
dl = -dl; 


