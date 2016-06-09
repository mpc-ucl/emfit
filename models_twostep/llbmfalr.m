function [l,dl] = llbmfalr(x,D,mu,nui,doprior);
%
% Fit only the model-free component of the model in Daw et al. 2011. Note this
% version assumes that betas at first and second levels take on the same value,
% and that the learning rates are the same, too. 
%
% Quentin Huys, 2015 
% www.quentinhuys.com/code.html 
% www.quentinhuys.com/pub.html
% qhuys@cantab.net

np = size(x,1);
if nargout==2; dodiff=1; else; dodiff=0;end


bmf  = exp(x(1));
al = 1./(1+exp(-x(2)));
la = 1./(1+exp(-x(3)));
rep = x(4)*eye(2);

Q1 = zeros(2,1);
Q2 = zeros(2,2);
dQ1da= zeros(2,1);
dQ2da= zeros(2,2);
dQeda= zeros(2,1);
dQ1dl = zeros(2,1);
dQedl = zeros(2,1);
drep = eye(2);

if doprior;
	l = -1/2*(x-mu)'*nui*(x-mu) - np/2*log(2*pi) - 1/2*log(1/det(nui));
	dl = - nui*(x-mu);
else;
	l=0;
	dl=zeros(np,1);
end

bb=20;
n=zeros(2);
for t=1:length(D.A);

	
	s=D.S(1,t); sp=D.S(2,t)-1;
	a=D.A(1,t); ap=D.A(2,t);
	r=D.R(1,t);

	if ~isnan(a) & ~isnan(ap); 

		if n(1,1)+n(2,2) > n(1,2)+n(2,1)
			Tr = .3+.4*eye(2);
		else
			Tr = .3+.4*(1-eye(2));
		end

		Qeff= bmf*Q1;
		if t>1 & exist('a1old'); Qeff= Qeff + rep(:,a1old); end

		lpa = Qeff;
		lpa = lpa - max(lpa);
		lpa = lpa - log(sum(exp(lpa)));
		l = l + lpa(a);
		pa = exp(lpa);

		lpap = bmf*Q2(:,sp);
		lpap = lpap - max(lpap);
		lpap = lpap - log(sum(exp(lpap)));
		l = l + lpap(ap);
		pap = exp(lpap);

		de1 = Q2(ap,sp)-Q1(a);
		de2 = r - Q2(ap,sp);


		if dodiff

			dl(1) = dl(1) + bmf*(Q1(a)-pa'*Q1 + Q2(ap,sp)-pap'*Q2(:,sp));

			dl(2) = dl(2) + bmf*(dQ1da(a) - pa'*dQ1da + dQ2da(ap,sp) - pap'*dQ2da(:,sp));

			dl(3) = dl(3) + dQedl(a)  - pa'*dQedl;

			% grad wrt al(1)
			dQ1da(a) = dQ1da(a) + al*(1-al)*(de1+la*de2) + al*(dQ2da(ap,sp)-dQ1da(a) + la*-dQ2da(ap,sp)); 
			dQ2da(ap,sp) = dQ2da(ap,sp) + al*(1-al)*de2 + al*-dQ2da(ap,sp);

			% grad wrt la
			dQ1dl(a) = dQ1dl(a) + al*(-dQ1dl(a) + la*(1-la)*de2);
			dQedl = bmf*dQ1dl;

			% grad wrt rep 
			if t>1 & exist('a1old');
				dl(4) = dl(4) + ((a==a1old) - pa'*drep(:,a1old));
			end

		end

		Q1(a)     = Q1(a)     + al*(de1 + la*de2);
		Q2(ap,sp) = Q2(ap,sp) + al*de2;

		n(sp,a) = n(sp,a)+1;
		a1old = a; 
	end

end
l=-l;
dl=-dl;
