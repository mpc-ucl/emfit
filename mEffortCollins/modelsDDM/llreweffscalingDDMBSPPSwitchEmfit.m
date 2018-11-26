function [l,dl, dsurr] = llreweffscalingDDMBSPPSwitchEmfit(x,D,mu,nui,doprior,options)

% This model consits of the analytical version of the drift
% diffusion model by Navarro % Fuss (2009). Its fits parameters for the starting point, 
% the boundary and non-decision time. The drift rate contains a standard model 
% for effort and reward evaluation, which assumes participants weigh both rewards and
% effort in their choices. The later model fits a weight for the rewards and the efforts. 

np = length(x);
dodiff= nargout==2;

spf = 1/(1+exp(-x(1)));  % parameter for starting point fraction
b = exp(x(2));           % parameter for boundary
betarew = exp(x(3));     % beta for reward
betaeff = exp(x(4));     % beta for effort
ndt = exp(x(5));         % parameter for non-decision time
pswitch = 1/(1+exp(-x(6)));

sp = spf*b;             % starting point is a fraction of the boundary

[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.rew; 
totalT = D.decisiontime;

effortCostLo = D.effortCostLo;%  set to - 0.2 in original task; 
effortCostHi = D.effortCostHi;%  set to - 1 in original task; 
effortCost = [effortCostLo effortCostHi]';


% if dogenerate
%     Vlow =  betarew*1    - betaeff*Effortcost2;     
%     for i = 1:5
%         reward = i+2; % 3,4,5,6,7
%         Vhigh(i) = betarew*reward- betaeff*Effortcost1;
%         % define drift rate as difference of value for high and low option
%         % "correct" boundary is lower boundary.
%         v = -1*(Vhigh(i)-Vlow);  
%         rew = 0; % not needed in this model
%         bscale = 0; % not needed in this model 
%         % make probability distributions with fitted parameters for
%         % possible time intervals to speed up simulations
%         [combined_t(i,:), combined_prob(i,:)]=make_cdfs_sp(v,b,sp,rew, bscale);      
%     end
%     %plotsprobs(combined_t, combined_prob) % to check if fitted
%     %distributions look reasonable
% end


for t=1:length(a)
	% define in terms of task 
	Vhigh = betarew*r(t)+betaeff*effortCostHi;
	Vlow =  betarew*1+betaeff*effortCostLo;
    % define drift rate as difference of value for high and low option
    % here defined such that drift rate is negative if goes to high value
    % option
    v = -1*(Vhigh-Vlow); 
    % only actually decision time is taking into account in DDM, thus
    % non-decision time is subtracted from recorded time 
    decT = totalT(t)-ndt; 
    [pt, dv, da, dz, dt] = wfpt_prep(b,v,sp,decT);
    
    

	if options.generatesurrogatedata==1
        rew = 0; 
        bscale = 0; 
        [combined_t(:), combined_prob(:)]=make_cdfs_sp(v,b,sp,rew, bscale); 
		[asurr(t), simTime(t)] = generateDataDDM(combined_t(:), combined_prob(:), ndt);
        if r(t) < 5 && asurr(t) == 2
            if rand<pswitch 
                asurr(t) = 1; 
            end
        end
    else
        if r(t) == 3 || r(t) == 4; 
            if a(t) == 1 % low
                l = l+log(pt(1)+pt(2)*pswitch);
            elseif a(t) == 2 % high
                l = l+log(pt(2)*(1-pswitch));
            end
        else 
            l = l+log(pt(a(t))); 
    end
    
        
    end
    
    
end
if options.generatesurrogatedata==1
	dsurr.a = asurr; 
    dsurr.simTime = simTime; 
end

dl = 0; 
l = -l; 
end
