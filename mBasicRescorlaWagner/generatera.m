function [a,r] = generatera(pa,t);

a = sum(rand>cumsum([0 pa]));

r = sign(sin(t/30*2*pi)) * (2*(rand<0.8)-1);
if a==2; r = - r;end

