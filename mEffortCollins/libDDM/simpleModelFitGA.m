options = optimoptions('ga', 'PlotFcn', @gaplotbestf,'FunctionTolerance', 10^(-5), 'MaxGenerations', 800, 'MaxStallGenerations', 200);
lb = [-6, -3, -3, -3, -3];
ub = [6, 3, 3, 3, 0];
for i = 1:Nsj
    D = Data(i); 
    save D; 
    [x,fval] = ga(@llreweffscalingDDMBSPGA,5,[],[],[],[],lb,ub,[],[],options);
    Data(i).estParams = x; 
    Data(i).lik = fval; 
end

 