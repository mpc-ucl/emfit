function [a,r] = genrw(x,T)

beta = exp(x(1));            % sensitivity to reward          
eps = 1./(1+exp(-x(2)));     % learning rate

pa=zeros(2,1);
Q=zeros(2,1); 

for t=1:T;
	l0 = max(Q);
	la = Q - l0 - log(sum(exp(Q-l0)));
	pa = exp(la);

	a(t) = (rand>pa(1))+1;

	r(t) = sign(sin(t/30*2*pi)) * (2*(rand<0.8)-1);
	if a(t)==2; r(t) = - r(t);end

	er = beta * r(t);

	Q(a(t)) = Q(a(t)) + eps * (er - Q(a(t)));  

end


