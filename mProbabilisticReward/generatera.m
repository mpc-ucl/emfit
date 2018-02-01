function [a,r] = generatera(pa,s,Z);

error('This isn''t properly implemented yet - and differs between different versions of the task, so be careful!');

	a = sum(rand>cumsum([0 pa']));

	if Z.I(a,s)
		r = rand<Z.prc(a);
	else
		r = 0;
	end
