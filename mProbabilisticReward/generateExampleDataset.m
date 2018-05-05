function Data=generateExampleDataset(Nsj,resultsDir)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for probabilistic reward task using the
% standard model llbelq0. 
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for probabilistic reward task\n')

options.generatesurrogatedata=1; 

T = 300; 
for sj=1:Nsj; 
	Data(sj).ID = sprintf('Subj %i',sj);
	
	Data(sj).a = zeros(T,1);					% preallocate space
	Data(sj).s = (rand(T,1)>0.5)+1;			% randomise stimuli 
	Data(sj).r = zeros(T,1);					% preallocate space

	Data(sj).bias = (rand>0.5)+1;				% rich stimulus 
	Data(sj).Nch = T; 							% length 
	Data(sj).I = eye(2); 						% instructions 
	Data(sj).prc = [0.78 0.33]';				% reward prob for rich/lean stimulus
	if Data(sj).bias==2; Data(sj).prc = 1-Data(sj).prc;end

	% realistic random parameters 
	Data(sj).trueParam = [1.3 0.3 -2 0.1 -0.1 ]'+.5*randn(5,1);

	% generate choices A, state transitions S and rewards R 
	[foo,foo,dsurr] = llbgelq0(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).a = dsurr.a;
	Data(sj).r = dsurr.r;
	Data(sj).trueModel='llbgelq0';

end

fprintf('Saved example dataset as Data.mat\n');
save([resultsDir filesep 'Data.mat','Data');
