function batchParameterPlots(Data,models,resultsDir,bestmodel)

nModls = length(models);
Nsj = length(Data);

nfig=0; 
eval(['mkdir ' resultsDir '/figs']);

%--------------------------------------------------------------------
% plot parameters of best model 
%--------------------------------------------------------------------

figure(2);clf;
for mdl=1:nModls
	R.(models(mdl).name) = load(sprintf('%s/%s',resultsDir,models(mdl).name));
end

for j=1:3

	if j==1
		% posterior parameters 
		par = R.(models(bestmodel).name).E;	
		yl = models(bestmodel).parnames_untr;
		yl{1} = {'EM-MAP',yl{1}};
	elseif j==2
		% transformed posterior parameters 
		R = modelConvertParams(R,models);	
		par = R.(models(bestmodel).name).parConvert;
		yl = models(bestmodel).parnames;
		yl{1} = {'EM-MAP',yl{1}};
	elseif j==3
		% lightly regularized ML parameters - closest to the data
		par = R.(models(bestmodel).name).stats.EMAP0; 
		npar= models(bestmodel).npar;
		yl = models(bestmodel).parnames_untr;
		yl{1} = {'MAP',yl{1}};
	end

	npar = size(par,1);
	m = mean(par,2);
	s = std(par,[],2)/sqrt(Nsj);

	for k=1:npar
		subplot(3,npar,k+(j-1)*npar)
		mybar(m(k),.7);
		hon
		plot(linspace(0.7,1.3,Nsj),par(k,:),'k.');
		mydeb(0,m(k),s(k));
		hof
		set(gca,'xtick',[]);
		ylabel(yl{k});
	end

end
myfig(gcf,[resultsDir filesep 'figs/Parameters']);

%--------------------------------------------------------------------
% try plotting parameters against true values if ran on generated data
%--------------------------------------------------------------------

if isfield(Data,'trueModel');

	try 
		parTrue = [Data.trueParam];
		parFitEM = R.(Data(1).trueModel).E;
		parFitML = R.(Data(1).trueModel).stats.EML;

		truemodelid = find(strcmp({models.name},Data(1).trueModel));

		npar = size(parTrue,1);

		figure(3);clf;
		for k=1:npar
			subplot(1,npar,k)
			plot(parTrue(k,:),parTrue(k,:),'r-','linewidth',2);
			hon
			h(1) = plot(parTrue(k,:),parFitEM(k,:),'b+','markersize',20);
			h(2) = plot(parTrue(k,:),parFitML(k,:),'k.','markersize',20);
			hof
			xlabel('True Param Value');
			if k==1; ylabel('inferred param value');end
			if k==1; legend(h,{'EM-MAP','ML'},'location','northwest');end
			title(models(truemodelid).parnames_untr(k));
			set(gca,'fontsize',14);
			mytightaxes
		end

		myfig(gcf,[resultsDir filesep 'figs/ParameterRecovery']);
	catch 
		fprintf('Plotting true against inferred parameters failed.\n');
	end
end

