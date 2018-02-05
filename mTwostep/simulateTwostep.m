function [a,spr] = simulateTwostep(pa,s,trans,rewprob)

if ~exist('rewprob'); 
	a = sum(rand>cumsum([0 pa']));
	if trans; if a==1; spr = 1; else spr = 2; end
	else    ; if a==1; spr = 2; else spr = 1; end
	end
else
	a = sum(rand>cumsum([0 pa']));
	spr = rewprob(s,a)>rand;
end
