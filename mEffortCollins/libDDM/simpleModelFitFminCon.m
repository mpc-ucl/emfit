% model estimates for individuals subjects

options = optimset('Display','off');
warning off all
lb = [-6, -3, -3, -3, -3];
ub = [6, 3, 3, 3, 0];

for sj=1:Nsj; 
	
	fstr=str2func('llreweffscalingDDMBSP');	% turn variable into function call 
	x0 = zeros(5,1);             	% initialisation for fminunc 
	mu = x0;						% prior mean 
	nu = eye(5)*10; 	% prior variance 
	doprior = 1; 					% add a prior for regularization or not
    llopt.generatesurrogatedata = 0;
    data = Data(sj); 
    [x,fval] = fmincon(@(x)fstr(x,data,mu,inv(nu),doprior,llopt),x0,[],[],[],[],lb,ub,[],options);    
	Data(sj).estParams = x;
	Data(sj).llik = fval;
	Data(sj).BIC = bayesianInformationCriterion(fval, 60, 5); 

end


