function Data=generateExampleDataset(Nsj)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for effort task using the
% standard model. 
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for effort task\n')

options.generatesurrogatedata=1; 

load TrialSeq.mat; 

T = 60; 
for sj=1:Nsj; 
	Data(sj).ID = sprintf('Subj %i',sj);
	Data(sj).Nch = T; 							% length 
	Data(sj).a = zeros(T,1);					% preallocate space
	Data(sj).rew = TrialSeq(:,3);				% high reward options 
	Data(sj).effortCostLo = - 0.2; 			% standard setting for 20 button presse
	Data(sj).effortCostHi = - 1; 				% standard setting for 100 button presse

	% realistic random parameters 
	Data(sj).trueParam = [0.8 2.3]'+.6*randn(2,1);

	% generate surrogate behavioural data 
	[foo,foo,dsurr] = llrewardeffort(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).a = dsurr.a;
	Data(sj).trueModel='llrewardeffort';

end

