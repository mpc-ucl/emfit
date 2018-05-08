function surrogateDataPlots(Data,models,SurrogateData,bestmodel,fitResults)

nModls = length(models);
Nsj = length(Data);

nfig=0; 

mkdir figs 

%--------------------------------------------------------------------
% compare with surrogate data 
%--------------------------------------------------------------------
% either generate new data: 
nfig=nfig+1; figure(nfig);clf;

cr = [1 1 2 2];

ns = zeros(max([Data.Nch]),4);
as = zeros(max([Data.Nch]),4);
bs = zeros(max([Data.Nch]),4);
nns = zeros(max([Data.Nch]),4,nModls);
for sj=1:Nsj
	a = Data(sj).a; 
	s = Data(sj).s; 
	for ss=1:4
		i = s==ss; 
		as(1:sum(i),ss,sj) = a(i)==1;
		ns(1:sum(i),ss) = ns(1:sum(i),ss)+1;

		pc(ss,sj) = sum(a(i)==cr(ss))/sum(i);
	end

	for mdl=1:nModls;
		a = [SurrogateData(sj).(models(mdl).name).a]; 
		nsample = numel(SurrogateData(sj).(models(mdl).name)); 
		a = reshape(a,size(a,2)/nsample,nsample);
		b = mean(a==1,2);
		for ss=1:4
			i = s==ss; 
			bs(1:sum(i),ss,mdl,sj) = b(i);
			nns(1:sum(i),ss,mdl) = nns(1:sum(i),ss,mdl)+1;

			pcs(ss,sj,mdl) = mean(sum(a(i,:)==cr(ss))/sum(i));
		end
	end

end
mas = sum(as,3)./ns;
mbs = sum(bs,4)./nns;

Ti = {'Go to win','Go to avoid','Nogo to win','Nogo to avoid'};
ssi = [1 3 2 4];

subplot(1,5,1)
	mybar(sum(pc(ssi,:)')/Nsj,.7);
	hon
	xx = [1:4]'*ones(1,nModls) + ones(4,1)*linspace(-.3,.3,nModls);
	plot(xx,sq(sum(pcs(ssi,:,:),2)/Nsj),'.-','markersize',15,'linewidth',1)
	hof
	xlim([.5 4.5]);
	ylabel('Probability correct');
	set(gca,'xticklabel',Ti(ssi),'xticklabelrotation',30);
		
for ss=1:4
	subplot(1,5,1+ss);
		plot(mas(:,ss),'k','linewidth',3);
		hon
		plot(sq(mbs(:,ss,:)),'linewidth',2)
		hof
		ylim([0 1]);
		title(Ti{ss});
		xlabel('Trial');
		if ssi(ss)==1; ylabel('Probability Go');end
end

le = {models.name}; 
le = {'Data',le{:}};
legend(le,'location','best'); 

myfig(gcf,sprintf('%s/figs/SurrogateDataPlots',fitResults));

