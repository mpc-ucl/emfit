function [a,r] = generatera(pa,s,prc,I);

if isnan(s)
	a = NaN;
	r = NaN;
else
	a = sum(rand>cumsum([0 pa']));
	if I(a,s)
		r = rand<prc(a);
	else
		r = 0;
	end
end
