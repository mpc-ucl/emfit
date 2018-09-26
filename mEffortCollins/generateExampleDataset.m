function Data=generateExampleDataset(Nsj,resultsDir,allResults)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for effort task by Gold, J.
% M.; Strauss, G. P.; Waltz, J.  A.; Robinson, B. M.; Brown, J. K. & Frank, M.
% J. Negative symptoms of schizophrenia are associated with abnormal effort-cost
% computations. Biol Psychiatry, 2013, 74, 130-136
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for effort task\n')

options.generatesurrogatedata=1; 

load TrialSeq.mat; 	% this is the fixed sequence of reward options presented in
							% the task. 

T = 60; 
for sj=1:Nsj; 
	Data(sj).ID = sprintf('Subj %i',sj);
	Data(sj).Nch = T; 					    % length 
	Data(sj).a = zeros(T,1);                % preallocate space
    Data(sj).decisionTime = zeros(T,1);     % preallocate space
	Data(sj).rew =allResults(sj).rewardTrace;		    % high reward options 
	Data(sj).effortCostLo = - 0.2; 			% standard setting for 20 button presse
	Data(sj).effortCostHi = - 1; 		    % standard setting for 100 button presse

	% realistic random parameters 
    
	%Data(sj).trueParam = [0.2 0.7 -0.5 0.8 -0.5 ]'+.6*randn(5,1);
    Data(sj).trueParam = allResults(sj).llreweffscalingDDMBSP.params;
    
    %Data(sj).trueParam = rand(5,1)*8-4';
	%  generate surrogate behavioural data 
	[foo,foo,dsurr] = llreweffscalingDDMBSP(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).a = dsurr.a;
    Data(sj).decisionTime = dsurr.simTime;  
	Data(sj).trueModel='llreweffscalingDDMBSP';

end

fprintf('Saved example dataset as Data.mat\n');
save([resultsDir filesep 'Data.mat'],'Data');
