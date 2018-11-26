
%Uncomment beneath if want to use model with or without pswitch and which
%optimization algorithm to use (Pswitch might only work with fminunc for
%now)
%Data=generateExampleDataset(Nsj,resultsDir, allResults); % for model without Pswitch
Data=generateExampleDatasetPSwitch(Nsj,resultsDir, allResults); % for model with Pswitch
simpleModelFit;             % uses fminunc 
%simpleModelFitFminCon;     % uses fmincon
%simpleModelFitGA;          % uses GA

paramsEst = [];
for i=1:Nsj
    paramsEst(i,:) = Data(i).estParams';
end

paramsTrue = [];
for i=1:Nsj
    paramsTrue(i,:) = Data(i).trueParam';
end

figure; 

subplot(2,3,1)
i = 1;
plot(paramsTrue(:,i), paramsEst(:,i),'.'); 
rho = corr(paramsTrue(:,i), paramsEst(:,i));
xlabel('true params'); 
ylabel('estimated params');
title(sprintf('starting point: rho = %.4f', rho)); 
xlim([-2 3])
ylim([-2 3])


subplot(2,3,2)
i = 2;
plot(paramsTrue(:,i), paramsEst(:,i),'.'); 
rho = corr(paramsTrue(:,i), paramsEst(:,i));
xlabel('true params'); 
ylabel('estimated params');
title(sprintf('boundary: rho = %.4f', rho)); 
xlim([-2 3])
ylim([-2 3])

subplot(2,3,3)
i = 5;
plot(paramsTrue(:,i), paramsEst(:,i),'.'); 
rho = corr(paramsTrue(:,i), paramsEst(:,i));
xlabel('true params'); 
ylabel('estimated params');
title(sprintf('non-decision time: rho = %.4f', rho)); 
xlim([-2 3])
ylim([-2 3])

subplot(2,3,4)
i = 3;
plot(paramsTrue(:,i), paramsEst(:,i),'.'); 
rho = corr(paramsTrue(:,i), paramsEst(:,i));
xlabel('true params'); 
ylabel('estimated params');
title(sprintf('betarew: rho = %.4f', rho)); 
xlim([-2 3])
ylim([-2 3])

subplot(2,3,5)
i = 4;
plot(paramsTrue(:,i), paramsEst(:,i),'.'); 
rho = corr(paramsTrue(:,i), paramsEst(:,i));
xlabel('true params'); 
ylabel('estimated params');
title(sprintf('betaeff: rho = %.4f', rho)); 
xlim([-2 3])
ylim([-2 3])

subplot(2,3,6)
i = 6;
plot(paramsTrue(:,i), paramsEst(:,i),'.'); 
rho = corr(paramsTrue(:,i), paramsEst(:,i));
xlabel('true params'); 
ylabel('estimated params');
title(sprintf('pswitch: rho = %.4f', rho)); 
xlim([-2 3])
ylim([-2 3])