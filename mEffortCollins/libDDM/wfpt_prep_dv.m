function pt = wfpt_prep_dv(b,v,dv,sp,time)
% Prepares data to calculate the derivatives for the drift rate parameters 
% of the probability for Wiener first passage time

pt = []; 
err = 10^(-29);


% Probability of making low choice (1)
% Drift rate is defined as negative if value for high options is better. 
% Thus, sign for drift rate must be swapped for low choice.
% Starting point is defined, such that 1-sp is the distance to boundaries
% for low choice
pt(1) = wfpt_dv(time,-v,-dv,b,b-sp,err); 


% Probability of making high choice (2)
pt(2) = wfpt_dv(time,v,dv,b,sp,err); 

end