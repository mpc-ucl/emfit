function pt = wfpt_prep(b,v,sp,time)
% Prepares data to calculate the probability for Wiener first passage time
% according to Navarro et al 2009.
pt = []; 
err = 10^(-29);

% Probability of making low choice (1)
% Drift rate is defined as negative if value for high options is better. 
% Thus, sign for drift rate must be swapped for low choice.
% Starting point is defined, such that 1-sp is the distance to boundaries
% for low choice
pt(1) = wfpt(time,-v,b,b-sp,err); 


% Probability of making high choice (2)
pt(2) = wfpt(time,v,b,sp,err); 

end         