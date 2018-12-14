function [l, dl, dsurr] = llreweffscalingDDMBScaledSPPSwitch(x,D,mu,nui,doprior,options)

% This model consits of the analytical version of the drift
% diffusion model by Navarro % Fuss (2009). Its fits parameters for the starting point, 
% the starting boundary, scaling boundary and non-decision time. The drift rate contains a standard model 
% for effort and reward evaluation, which assumes participants weigh both rewards and
% effort in their choices. The later model fits a weight for the rewards and the efforts. 

np = length(x);
dodiff= nargout==2;

spf = 1/(1+exp(-x(1)));  % parameter for starting point fraction
sb = exp(x(2));          % parameter for starting boundary
betarew = exp(x(3));     % beta for reward
betaeff = exp(x(4));     % beta for effort
ndtime = 1/(1+exp(-x(5)));         % parameter for non-decision time
btrialscale = 1/(1+exp(-x(6))); % parameter for scaling the boundary
pswitch = 1/(1+exp(-x(7)));

ndt = ndtime*0.7;

btrialscaled = btrialscale * (sb/60); % to make sure that boundary cannot get negative 

[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.rew; 
totalT = D.decisiontime;

effortCostLo = D.effortCostLo;%  set to - 0.2 in original task; 
effortCostHi = D.effortCostHi;%  set to - 1 in original task; 
effortCost = [effortCostLo effortCostHi]';





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
%     if decT < 0
%         decT = 10e-2;
%     end
    b = sb-btrialscaled*t; 
    sp = spf*b;              % starting point is a fraction of the boundary
    [pt, dv, da, dz, dt, dv_ps, dz_ps, da_ps, dt_ps]= wfpt_prep(b,v,sp,decT);
  
	if options.generatesurrogatedata==1
        rew = 0; 
        bscale = 0; 
        [combined_t(:), combined_prob(:)]=make_prob_dist(v,b,sp,rew,bscale);      
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
    
    if dodiff
       % derivative of starting point
       if r(t) == 3 || r(t) == 4; 
            if a(t) == 1 % low
                dl(1) = dl(1)+(1/(pt(1)+pt(2)*pswitch))*(dz_ps(1)*(-1)*(spf*(1-spf))*b+dz_ps(2)*(spf*(1-spf))*b*pswitch);
            elseif a(t) == 2 % high
                dl(1) = dl(1)+dz(2)*(spf*(1-spf))*b;
            end
       else 
           if a(t) == 1
               dl(1) = dl(1)+dz(a(t))*(-1)*(spf*(1-spf))*b; 
           else
               dl(1) = dl(1)+dz(a(t))*(spf*(1-spf))*b;
           end
       end
       
       % derivative of boundary
       if r(t) == 3 || r(t) == 4; 
            if a(t) == 1 % low
                dl(2) = dl(2)+(1/(pt(1)+pt(2)*pswitch))*((da_ps(1)*b+dz_ps(1)*(b-spf*b)+(da_ps(2)*b+dz_ps(2)*spf*b)*pswitch));
            elseif a(t) == 2 % high
                dl(2) = dl(2)+(da(a(t))*b+dz(a(t))*spf*b);
            end
       else 
           if a(t) == 1
               dl(2) = dl(2)+(da(a(t))*b+dz(a(t))*(b-spf*b));
           else 
               dl(2) = dl(2)+(da(a(t))*b+dz(a(t))*spf*b);
           end
       end
           
       % derivative of betarew
        dvbr = -1*(betarew*r(t)-betarew*1); 
       if r(t) == 3 || r(t) == 4; 
           if a(t) == 1 % low
                dl(3) = dl(3)+(1/(pt(1)+pt(2)*pswitch))*(dv_ps(1)*-dvbr+dv_ps(2)*dvbr*pswitch);
            elseif a(t) == 2 % high
                dl(3) = dl(3)+dv(a(t))*dvbr; 
           end     
       else    
           if a(t) == 1
               dl(3) = dl(3)+dv(a(t))*-dvbr;  
           elseif a(t) == 2
               dl(3) = dl(3)+dv(a(t))*dvbr; 
           end
       end
       
       % derivative of betaeff
       dvbe = -1*(betaeff*effortCostHi-betaeff*effortCostLo); 
       if r(t) == 3 || r(t) == 4; 
           if a(t) == 1 % low
                dl(4) = dl(4)+(1/(pt(1)+pt(2)*pswitch))*(dv_ps(1)*-dvbe+dv_ps(2)*dvbe*pswitch);
            elseif a(t) == 2 % high
                dl(4) = dl(4)+dv(a(t))*dvbe; 
           end     
       else    
           if a(t) == 1
               dl(4) = dl(4)+dv(a(t))*-dvbe;  
           elseif a(t) == 2
               dl(4) = dl(4)+dv(a(t))*dvbe; 
           end
       end
       
       % derivative of non-decision time
       if r(t) == 3 || r(t) == 4; 
           if a(t) == 1 % low
                dl(5) = dl(5)+(1/(pt(1)+pt(2)*pswitch))*(dt_ps(1)*(-0.7*(ndtime*(1-ndtime)))+dt_ps(2)*(-0.7*(ndtime*(1-ndtime)))*pswitch);
            elseif a(t) == 2 % high
                dl(5) = dl(5)+dt(a(t))*(-0.7*(ndtime*(1-ndtime))); 
           end     
       else   
           dl(5) = dl(5)+dt(a(t))*(-0.7*(ndtime*(1-ndtime))); 
       end
       
        % derivative of boundary scaling
       dbts = (-1)*btrialscale*(1-btrialscale)*(sb/60)*t;
       if r(t) == 3 || r(t) == 4; 
            if a(t) == 1 % low
                dl(6) = dl(6)+(1/(pt(1)+pt(2)*pswitch))*((da_ps(1)*dbts+dz_ps(1)*(dbts-spf*dbts)+(da_ps(2)*dbts+dz_ps(2)*spf*dbts)*pswitch));
            elseif a(t) == 2 % high
                dl(6) = dl(6)+(da(a(t))*b+dz(a(t))*spf*dbts);
            end
       else 
           if a(t) == 1
               dl(6) = dl(6)+(da(a(t))*dbts+dz(a(t))*(dbts-spf*dbts));
           else 
               dl(6) = dl(6)+(da(a(t))*dbts+dz(a(t))*spf*dbts);
           end
       end
       
       % derivative of pswitch
       if r(t) == 3 || r(t) == 4; 
            if a(t) == 1; 
                dl(7) = dl(7)+(1/(pt(1)+pt(2)*pswitch))*pt(2)*pswitch*(1-pswitch);         
            elseif a(t) == 2 % high
                dl(7) = dl(7)+(1/(1-pswitch))*(-1)*pswitch*(1-pswitch);
            end
       else
           dl(7) = dl(7)+0; 
       end       
    end
    
end

l = -l; 
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = asurr; 
    dsurr.simTime = simTime; 
end
end


