function [asurr, simTime] = simulateEffortDDM(combined_t, combined_prob,ndt)

dt = 0.01;

% make cdf of probability distribution over low and high choices
cdf_prob = cumsum(combined_prob); 

% get time corresponding to randomly choosen prob value of cdf
r1 = rand(1);
A = cdf_prob < r1;   
    try 
        simTime = abs(combined_t(sum(A)+1)-dt/2)+ndt;
    catch 
        % in case highest interval was choosen
        simTime = abs(combined_t(sum(A))-dt/2)+ndt;
    end
        
% get choice corresponding to same randomly choosen prob value of cdf         
    try 
        if combined_t(sum(A)) > 0 % if time is positive, choice is high 
            asurr = 2;
        else
            asurr = 1;
        end
    catch
        % if time is 0
        asurr = 2; 
    end
    
end
   