function [a,r] = generatera(pa,s);

p = 0.8; 
c = 0.4;
ago = 1; 

a = sum(rand>cumsum([0 pa]));

if s==1	% go to win 
	if a == ago; 	r = c*(rand<p);
	else 			r = 0; 
	end
elseif s==2; 	% go to avoid 
	if a == ago; 	r = 0; 
	else 			r = -c*(rand<p);
	end
elseif s==3		% nogo to win 
	if a == ago; 	r = 0; 
	else 			r = c*(rand<p);
	end
elseif s==4; 	% nogo to avoid 
	if a == ago; 	r = -c*(rand<p);
	else 			r = 0; 
	end
end
