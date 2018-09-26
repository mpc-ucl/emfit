
paramRange1 = -6:0.1:6;
paramRange2 = -6:0.1:6;

%D = Data(2); 
x = [0.2000    0.7000   -0.5000    0.8000   -0.5000];
pn1 = 1;
pn2 = 2; 
param1 = 'starting point'; 
param2 = 'boundary';

for i = 1:length(paramRange1)
    x(pn1) = paramRange1(i); 
    for j = 1:length(paramRange2)
        x(pn2) = paramRange2(j); 

        mu = 0; 
        nui = 0; 
        doprior = 0; 

        dodiff =1; 

        np = length(x);


        spf = 1/(1+exp(-x(1)));  % parameter for starting point fraction
        b = exp(x(2));           % parameter for boundary
        betarew = exp(x(3));     % beta for reward
        betaeff = exp(x(4));     % beta for effort
        ndt = exp(x(5));         % parameter for non-decision time


        sp = spf*b;              % starting point is a fraction of the boundary

        [l,dl] = logGaussianPrior(x,mu,nui,doprior);

        a = D.a;%D.asurr(:,1)+1; 
        r = D.rew; 
        totalT = D.decisionTime;%D.timesurr(:,1);

        effortCostLo = D.effortCostLo;%  set to - 0.2 in original task; 
        effortCostHi = D.effortCostHi;%  set to - 1 in original task; 
        effortCost = [effortCostLo effortCostHi]';



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
            l = l+log(pt(a(t)));
        end

            if dodiff
               % derivative of starting point
               ptdz = wfpt_prep_dz(b,v,sp,decT);
               dl(1) = ptdz(a(t)); 
               % derivative of boundary
                ptda = wfpt_prep_da(b,v,sp,decT);
                dl(2) = ptda(a(t)); 
               % derivative of betarew
               dvbr = -1*(betarew*r(t)-betarew*1); 
               pt = wfpt_prep_dv(b,v,dvbr,sp,decT);
               dl(3) = pt(a(t)); 
               % derivative of betaeff
               dvbe = -1*(betaeff*effortCostHi-betaeff*effortCostLo); 
               pt = wfpt_prep_dv(b,v,dvbe,sp,decT);
               dl(4)= pt(a(t)); 
               % derivative of non-decision time
            end
        l = -l; 
        dl = -dl;
        dlpn1(i,j) = dl(pn1); 
        dlpn2(i,j) = dl(pn2); 
        lik2(i,j) = l;
        end

       
end



figure; h = surf(paramRange1, paramRange2, lik2, 'FaceAlpha',1);set(h,'LineStyle','none');
xlabel(param2)
ylabel(param1)
zlabel('neg likelihood')

figure; h = surf(paramRange1, paramRange2, dlpn1, 'FaceAlpha',1);set(h,'LineStyle','none');
xlabel(param2)
ylabel(param1)
zlabel('der starting point')

figure; h = surf(paramRange1, paramRange2, dlpn2, 'FaceAlpha',1);set(h,'LineStyle','none');
xlabel(param2)
ylabel(param1)
zlabel('der boundary')

figure; plot(paramRange1, lik2(60,:));
xlabel(param2)
ylabel('neg likelihood')

figure; plot(paramRange1, dlpn1(55,:));
xlabel(param2)
ylabel('dersp')

figure; plot(paramRange1, dlpn2(:,100));
xlabel(param1)
ylabel('der boundary')

figure; plot(paramRange1, dlpn2(60,:));
xlabel(param2)
ylabel('der boundary')
