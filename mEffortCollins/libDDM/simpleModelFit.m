% model estimates for individuals subjects

options = optimoptions('fminunc','Algorithm','trust-region','CheckGradients',false,'Display','none');

for sj=1:Nsj; 
	
	fstr=str2func('llreweffscalingDDMBSPPSwitchEmfit');	% turn variable into function call 
    %fstr=str2func('llreweffscalingDDMBSPPSwitchEmfit');
	x0 = zeros(6,1);             	% initialisation for fminunc 
	mu = x0;						% prior mean 
	nu = eye(6)*10; 	% prior variance 
	doprior = 1; 					% add a prior for regularization or not
    llopt.generatesurrogatedata = 0;
    data = Data(sj); 
    sto = 1; 
    while sto == 1
        try
        [x,fval,exitflag] = fminunc(@(x)fstr(x,data,mu,inv(nu),doprior,llopt),x0,options);
        sto = 0; 
        catch 
            x0 = zeros(6,1)+.6*randn(6,1);
        end
    end
    
	Data(sj).estParams = x;
	Data(sj).llik = fval;
    Data(sj).exitflag = exitflag; 
	Data(sj).BIC = bayesianInformationCriterion(fval, 60, 6); 

end

