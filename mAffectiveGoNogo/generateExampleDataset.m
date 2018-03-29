function Data=generateExampleDataset(Nsj)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for affective Go/Nogo task using the
% standard model llbaepqx.m
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for affective Go/Nogo task\n')

options.generatesurrogatedata=1; 

T = 160; 
for sj=1:Nsj; 
	Data(sj).ID = sprintf('Subj %i',sj);
	
	Data(sj).a = zeros(1,T);					% preallocate space
	Data(sj).r = zeros(1,T);					% preallocate space
	Data(sj).s = randi(4,T,1);					% randomise stimuli 

	Data(sj).Nch = T; 							% length 

	% realistic random parameters 
	Data(sj).trueParam = [3.4 -2.5 -1 1.3 4.2 ]'+randn(5,1);

	% generate choices A, state transitions S and rewards R 
	[foo,foo,dsurr] = llbaepqx(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).a = dsurr.a;
	Data(sj).r = dsurr.r;
	Data(sj).trueModel='llbaepqx';

end

fprintf('Saved example dataset as Data.mat\n');
save Data.mat Data; 



