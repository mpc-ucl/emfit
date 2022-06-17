function [a,r] = generatera(pa,wordval,avatval);

a = (rand>pa(1))+1;

if ( a == 1 & (wordval ==  avatval) ) |  ...
	( a == 2 & (wordval ~=  avatval) );
	  r = 1; 
else r = -1; 
end

