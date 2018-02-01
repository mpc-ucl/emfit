function [l,dl,dsurr] = llrw(x,D,mu,nui,doprior,options);
% 
% function [l,dl,a,r] = llrw(x,D,mu,nui,doprior,options);
% 
% Likelihood function for simple Rescorla-Wagner Q learning model. It computes
% the probabilities of actions 
%
%       p(D.A(t)) = exp(beta*Q_t(D.A(t),D.S(t))) / sum_a exp(beta*Q_t(a,D.S(t)))
%
% with 
%
%       Q_t(D.A(t),D.S(t)) = Q_{t-1}(D.A(t),D.S(t)) +  epsilon*PE_t
%       PE_t =  D.R(t) - Q_{t-1}(D.A(t),D.S(t))
%   
% It outputs the total log likelihood L of all the choices, in addition to the
% gradient wrt. to the paramters X. Note that the parameters X are transformed
% to lie on the real line.  It applies a Gaussian prior with mean MUSJ and
% inverse variance nui if DOPRIOR=1. 
% 
% Copyright Quentin Huys, 2018 www.quentinhuys.com qhuys@cantab.net

dodiff=nargout==2;
np = length(x);
beta 		= exp(x(1));				% sensitivity to reward          
alfa 		= 1./(1+exp(-x(2)));		% learning rate

if ~exist('options'); options.generatesurrogatedata=0;end

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

% initialize some variables 
pa=zeros(2,1);
Q=zeros(2,1); 
dQdb = zeros(2,1);
dQde = zeros(2,1);

if options.generatesurrogatedata==1
	dodiff=0;	% don't compute gradients 
	a = zeros(size(D.a));
	r = zeros(size(D.a));
else 	% extract data needed here from D 
	a = D.a; 
	r = D.r; 
end

% loop over all choices 
for t=1:length(a)

	% compute choice likelihood given parameter - softmax of Q 
	l0 = max(Q);
	la = Q - l0 - log(sum(exp(Q-l0)));
	pa = exp(la);

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pa',t);	 	% generate an action and reward 
	else
		l = l + la(a(t));								% accumulate likelihood 
	end

	er = beta * r(t);

	% compute gradients of likelihood wrt parameter 
	if dodiff
		dl(1) = dl(1) + dQdb(a(t)) - pa'*dQdb;
		dl(2) = dl(2) + dQde(a(t)) - pa'*dQde;

		dQdb(a(t)) = (1-alfa)*dQdb(a(t)) + alfa*er;
		dQde(a(t)) = (1-alfa)*dQde(a(t)) + alfa*(1-alfa)*(er-Q(a(t)));
	end

	% update Q values 
	Q(a(t)) = Q(a(t)) + alfa * (er - Q(a(t)));  

end

% use f*MIN*unc, so minimize negative log likelihood 
l  = -l ;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

return 
