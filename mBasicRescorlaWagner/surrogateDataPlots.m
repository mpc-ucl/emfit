function surrogateDataPlots(Data,models,SurrogateData,bestmodel,fitResults)

nModls = length(models);
Nsj = length(Data);

nfig=get(gcf,'Number');

mkdir figs 

%--------------------------------------------------------------------
% compare with surrogate data 
%--------------------------------------------------------------------

fprintf('No surrogate data plots defined for basic RW model\n');
