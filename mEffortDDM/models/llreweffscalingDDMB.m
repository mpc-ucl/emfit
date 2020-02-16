function [l, dl, dsurr] = llreweffscalingDDMB(x,D,mu,nui,doprior,options)

% This model conists of the analytical version of the drift
% diffusion model by Navarro % Fuss (2009). It fits parameters for the the boundary and non-decision time. 
% The drift rate contains a standard model 
% for effort and reward evaluation, which assumes participants weigh both rewards and
% effort in their choices. The later model fits a weight for the rewards and the efforts. 

np = length(x);
dodiff= nargout==2;

spf = 0.5;               % starting point is fixed in the middle of the high and the low boundary
b = exp(x(1));           % parameter for boundary
betarew = exp(x(2));     % beta for reward
betaeff = exp(x(3));     % beta for effort
ndtime = 1/(1+exp(-x(4)));         % parameter for non-decision time

ndt = ndtime*0.7;        % ndt cannot be bigger than 0.7
sp = spf*b;              % starting point is a fraction of the boundary

[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.rew; 
totalT = D.decisiontime;

effortCostLo = D.effortCostLo;%  set to - 0.2 in original task; 
effortCostHi = D.effortCostHi;%  set to - 1 in original task; 
effortCost = [effortCostLo effortCostHi]';

if options.generatesurrogatedata==1
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
    [pt, dv, da, dz, dt] = wfpt_prep(b,v,sp,decT);
    
	if options.generatesurrogatedata==1
        rewidx = r(t)-2; % reward index 
        sto = 1;
        while sto == 1        
            [asurr(t), simTime(t)] = generateDataDDM(combined_t(rewidx,:), combined_prob(rewidx,:), ndt);
            if simTime(t) > 0.7 % fitting of ndt only works if no trials are slower than 0.7
                sto = 0;
            end
        end 
    else 
        l = l+log(pt(a(t)));
    end
    
    if dodiff
       % derivative of boundary
       if a(t) == 1
           dl(1) = dl(1)+(da(a(t))*b+dz(a(t))*(b-spf*b));
       else 
           dl(1) = dl(1)+(da(a(t))*b+dz(a(t))*spf*b);
       end
       % derivative of betarew
       dvbr = -1*(betarew*r(t)-betarew*1); 
       if a(t) == 1
           dl(2) = dl(2)+dv(a(t))*-dvbr;  
       elseif a(t) == 2
           dl(2) = dl(2)+dv(a(t))*dvbr; 
       end
       % derivative of betaeff
       dvbe = -1*(betaeff*effortCostHi-betaeff*effortCostLo); 
       if a(t) == 1
           dl(3) = dl(3)+dv(a(t))*-dvbe;
       elseif a(t) == 2
           dl(3) = dl(3)+dv(a(t))*dvbe;
       end
       % derivative of non-decision time
       dl(4) = dl(4)+dt(a(t))*(-0.7*(ndtime*(1-ndtime))); 
    end
    
end

l = -l; 
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = asurr; 
    dsurr.simTime = simTime; 
end
end


