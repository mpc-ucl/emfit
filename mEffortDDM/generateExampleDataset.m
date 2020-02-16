function Data=generateExampleDataset(Nsj,resultsDir)
% 
% Data = generateExampleDataset(Nsj)
% Berwian, I.M., Wenzel, J.G, Collins, A.G.E, Seifritz, E., Stephan, K.E., 
% Walter, H. & Huys, Q.J.M. (2020). Computational Mechanisms of Effort and 
% Reward Decisions in Patients With Depression and Their Association With 
% Relapse After Antidepressant Discontinuation. JAMA Psychiatry. 
% doi:10.1001/jamapsychiatry.2019.4971 
 
% Isabel Berwian & Quentin Huys 2020 www.quentinhuys.com 


fprintf('Generating example dataset for effort task\n')

options.generatesurrogatedata=1; 

load TrialSeq.mat; 	% this is a fixed sequence of reward options presented in
							% the task. 

T = 60; 
for sj=1:Nsj; 
	Data(sj).ID = sprintf('Subj %i',sj);
	Data(sj).Nch = T; 					    % length 
	Data(sj).a = zeros(T,1);                % preallocate space
    Data(sj).decisiontime = zeros(T,1);     % preallocate space
    Data(sj).rew = TrialSeq(:,3);   	    % high reward options 
	Data(sj).effortCostLo = - 0.2; 			% standard setting for 20 button presse
	Data(sj).effortCostHi = - 1; 		    % standard setting for 100 button presse

	% realistic random parameters 
	Data(sj).trueParam = [-0.15 0.9 -0.6 0.6 2.8 ]'+.3*randn(5,1);
  
	%  generate surrogate behavioural data 
	[foo,foo,dsurr] = llreweffscalingDDMBSP(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).a = dsurr.a;
    Data(sj).decisiontime = dsurr.simTime;  
	Data(sj).trueModel='llreweffscalingDDMBSP';

end

fprintf('Saved example dataset as Data.mat\n');
save([resultsDir filesep 'Data.mat'],'Data');
