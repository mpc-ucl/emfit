function pt = wfpt_prep_dz(b,v,sp,time)
% Prepares data to calculate the derivatives for the drift rate parameters 
% of the probability for Wiener first passage time

pt = []; 
err = 10^(-29);


% Probability of making low choice (1)
% Drift rate is defined as negative if value for high options is better. 
% Thus, sign for drift rate must be swapped for low choice.
% Starting point is defined, such that 1-sp is the distance to boundaries
% for low choice
pt(1) = wfpt_dz_low_choice(time,-v,b,sp,err); % 1-sp is taking inside the
% function, thus feed sp here


% Probability of making high choice (2)
pt(2) = wfpt_dz_high_choice(time,v,b,sp,err); 

end