function Data=generateExampleDataset(Nsj)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for pruning task using the
% full model llsrho2p.m 
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for Pruning task\n')

options.generatesurrogatedata=1; 

dmax=6; 
T = 100; 

Z.dmax = dmax; 
Z = precomputeParams(Z);
for sj=1:Nsj; 

	clear d dd ddd; 
	for k=1:T
		d  = randi(4)+2; 		% random depths from 3-6 
		dd.depths(k,1)  = d; 
		dd.choices(k,1:d) = ones(1,d);
		dd.states(k,1:d+1)  = repmat(randi(6),1,d+1);
	end
	dd.rewards = zeros(T,dmax);		% preallocate space

	ddd = extractValidTrials(dd);

	% realistic random parameters 
	trueParam(:,sj) = [-1.6 0 0.4 5 2 -0.3 ]'+randn(6,1);

	% generate choices A, state transitions S and rewards R 
	[foo,foo,dsurr] = llsrho2p(trueParam(:,sj),ddd,0,0,0,options); 

	an=dsurr.an;	
	sn=ddd.sn; 
	for k=1:T
		d = ddd.depths(k);
		D(sj).choices(k,1:d)   = Z.A(1:d,an(k),d-2);
		D(sj).states(k,1:d+1)  = Z.S(1:d+1,an(k),d-2,sn(k));
		D(sj).rewards(k,1:d)   = Z.R(1:d,an(k),d-2,sn(k));
	end

end

Data = extractValidTrials(D);
for sj=1:Nsj
	Data(sj).trueParam = trueParam(:,sj);
	Data(sj).trueModel='llsrho2p';
end

fprintf('Saved example dataset as Data.mat\n');
save Data.mat Data; 
