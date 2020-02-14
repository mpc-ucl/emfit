function [combined_t, combined_prob]=make_prob_dist(v,b,sp,rew,bscale) 
% This functions calculates probability distributions with fitted parameters 
% for possible time intervals to speed up simulation.

% define allowed error and time intervals  
err = 10^(-29); % error
dt = 0.01; % time intervals for simulation
t = [10^-80 dt:dt:7]; 
pl = zeros(length(t),1);
pu = zeros(length(t),1); 

% Probability of making low choice (pl)
% Drift rate is defined as negative if value for high options is better. 
% Thus, sign for drift rate must be swapped if low option is better. 
for i = 1:length(t)
    b = bscale*rew+b; % Boundary can vary depending on reward        
    pl(i) = wfpt(t(i),-v,b,b-sp,err); 
end

% Probability of making high choice (pu)
for i = 1:length(t)
    b = bscale*rew+b;
    pu(i) = wfpt(t(i),v,b,sp,err); 
end

% Uncomment this if want to check if probability distributions of low and
% high choices add up to one. 
%     total_prob_pl = sum(pl*dt)
%     total_prob_pu = sum(pu*dt)
%     total_prob = total_prob_pl+total_prob_pu
    
% get correct probability values    
pl_prob = pl*dt; 
pu_prob = pu*dt; 

% flipping time and probability for low choices on y axis to make later
% cdfs over both distributions at once.
flipped_prob_pl = fliplr(pl_prob'); 
combined_prob = [flipped_prob_pl'; pu_prob];
tleft = fliplr(-t); 
combined_t = [tleft, t];
 
end