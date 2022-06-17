function Data=generateExampleDataset(Nsj,resultsDir)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for affective Go/Nogo task using the
% standard model llbaepqx.m
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for self/other judgement task k\n')

options.generatesurrogatedata=1; 

wordval = [ -1 1 -1 1 -1 1 -1 1 -1 -1 1 -1 1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 1 -1 1 -1 -1 1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 1 -1 -1 1]';
avatval = [ 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 -1 1 -1 1 1 1 1 -1 1 1 1 -1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 1 1 -1 1 1 1 1 -1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 -1 1 1 1 1]';
avatid  = [1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 4 4 4 4 4 4 4 4 9 9 9 9 9 9 9 9 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 10 10 10 10 10 10 10 10 12 12 12 12 12 12 12 12 11 11 11 11 11 11 11 11]'; 

T = 96; 
for sj=1:Nsj; 
	Data(sj).ID = sprintf('Subj %i',sj);
	
	Data(sj).Nch = T; 							% length 
	Data(sj).a = zeros(T,1);					% preallocate space
	Data(sj).r = zeros(T,1);					% preallocate space
	Data(sj).wordval = wordval; 				% word values 
	Data(sj).avatval = avatval; 				% avatar values 
	Data(sj).avatid  = avatid ; 				% avatar id 

	% realistic random parameters 
	Data(sj).trueParam = [4.2 -5.3 0.8 0.3 ]'+randn(4,1);

	% generate choices A, state transitions S and rewards R 
	[foo,foo,dsurr] = lld6avb(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).a = dsurr.a;
	Data(sj).r = dsurr.r;
	Data(sj).trueModel='lld6avb';

end

fprintf('Saved example dataset as Data.mat\n');
save([resultsDir filesep 'Data.mat'],'Data');
