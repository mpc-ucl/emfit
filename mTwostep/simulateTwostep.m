function [a,spr] = simulateTwostep(pa,s,trans,rewprob)

a = sum(rand>cumsum(pa'))+1;

if nargin==4% try level 2 with rewprob 
	spr = rewprob(s,a)>rand;
else % if level 1 - when rewprob not provided 
	if trans; if a==1; spr = 1; else spr = 2; end
	else    ; if a==1; spr = 2; else spr = 1; end
	end
end
