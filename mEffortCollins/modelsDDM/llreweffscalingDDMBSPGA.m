function [l] = llreweffscalingDDMBSPGA(x)

% This model consits of the analytical version of the drift
% diffusion model by Navarro % Fuss (2009). Its fits parameters for the starting point, 
% the boundary and non-decision time. The drift rate contains a standard model 
% for effort and reward evaluation, which assumes participants weigh both rewards and
% effort in their choices. The later model fits a weight for the rewards and the efforts. 


generatesurrogatedata = 0; 
doprior = 0; 
mu = 0; 
nui = 1; 
load D


np = length(x);
dodiff= nargout==2;

spf = 1/(1+exp(-x(1)));  % parameter for starting point fraction
b = exp(x(2));           % parameter for boundary
betarew = exp(x(3));     % beta for reward
betaeff = exp(x(4));     % beta for effort
ndt = exp(x(5));         % parameter for non-decision time


sp = spf*b;              % starting point is a fraction of the boundary

[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.rew; 
totalT = D.decisionTime;

effortCostLo = D.effortCostLo;%  set to - 0.2 in original task; 
effortCostHi = D.effortCostHi;%  set to - 1 in original task; 
effortCost = [effortCostLo effortCostHi]';


if generatesurrogatedata==1
    dodiff=0;
    Vlow =  betarew*1+betaeff*effortCostLo;
    for i = 1:5
        Vhigh(i) = betarew*(i+2)+betaeff*effortCostHi;
        % define drift rate as difference of value for high and low option
        % "correct" boundary is lower boundary, therefore drift rate defined here as neg.  
        v = -1*(Vhigh(i)-Vlow); 
        rew = 0; % not needed in this model
        bscale = 0; % not needed in this model 
        % make probability distributions with fitted parameters for
        % possible time intervals to speed up simulations
        [combined_t(i,:), combined_prob(i,:)]=make_prob_dist(v,b,sp,rew,bscale);      
    end
    %plotsprobs(combined_t, combined_prob) % to check if fitted
    %distributions look reasonable
end


for t=1:length(a)
	% define in terms of task 
	Vhigh = betarew*r(t)+betaeff*effortCostHi;
	Vlow =  betarew*1+betaeff*effortCostLo;
    % define drift rate as difference of value for high and low option
    v = -1*(Vhigh-Vlow); 
    % only actually decision time is taking into account in DDM, thus
    % non-decision time is subtracted from recorded time 
    decT = totalT(t)-ndt; 
    pt = wfpt_prep(b,v,sp,decT);
    
    
	if generatesurrogatedata==1
        rewidx = r(t)-2; % reward index 
		[asurr(t), simTime(t)] = generateDataDDM(combined_t(rewidx,:), combined_prob(rewidx,:), ndt);
    else 
        l = l+log(pt(a(t)));
    end
    
    if dodiff
       % derivative of starting point
       % derivative of boundary
       % derivative of betarew
       dvbr = -1*(betarew*r(t)-betarew*1); 
       pt = wfpt_prep_dv(b,v,dvbr,sp,time);
       dl(3) = pt(a(t)); 
       % derivative of betaeff
       dvbe = -1*(betaeff*effortCostHi-betaeff*effortCostLo); 
       pt = wfpt_prep_dv(b,v,dvbe,sp,time);
       dl(4)= pt(a(t)); 
       % derivative of non-decision time
    end
    
end

l = -l; 


if generatesurrogatedata==1
	dsurr.a = asurr; 
    dsurr.simTime = simTime; 
end
end


