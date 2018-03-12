function Data=generateExampleDataset(Nsj)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for twostep task using the
% standard model llm2b2ar. 
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for twostep task\n')

load 'DawEtAl11RandomWalks.mat'

options.generatesurrogatedata=1; 

T = 201; 
for sj=1:Nsj; 
	Data(sj).ID = sprintf('Subj %i',sj);
	
	Data(sj).A = zeros(2,T);					% preallocate space
	Data(sj).S = [ones(1,T);zeros(1,T)];	% preallocate space
	Data(sj).R = zeros(1,T);					% preallocate space

	Data(sj).Nch = T; 							% length 

	Data(sj).trans = rand(1,T)<2/3;			% frequent or rare transition, predetermined
	Data(sj).rewprob = dawrandomwalks; 		% original reward random walk 

	% realistic random parameters 
	Data(sj).trueParam = [1.5 0.9 -0.2 -0.2 0.3 -0.3 0.2]'+.5*randn(7,1);

	% generate choices A, state transitions S and rewards R 
	[foo,foo,dsurr] = llm2b2alr(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).A = dsurr.A;
	Data(sj).S = dsurr.S;
	Data(sj).R = dsurr.R;
	Data(sj).trueModel='llm2b2alr';

end

fprintf('Saved example dataset as Data.mat');
save Data.mat Data; 




