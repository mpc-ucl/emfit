function surrogateDataPlots(Data,models,SurrogateData,bestmodel)

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

ns = zeros(40,4);
nns = zeros(40,4,nModls);
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

myfig(gcf,'figs/SurrogateDataPlots');

%--------------------------------------------------------------------
% plot parameters of best model 
%--------------------------------------------------------------------

nfig=nfig+1; figure(nfig);clf;
for mdl=1:nModls
	R.(models(mdl).name) = load(sprintf('fitResults/%s',models(mdl).name));
end

%try 

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
	myfig(gcf,'figs/Parameters');

%catch 
%	fprintf('No fits for model %s found, no parameters plotted \n',models(bestmodel).name);
%end



