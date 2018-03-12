function surrogateDataPlots(Data,models,SurrogateData,bestmodel)

nModls = length(models);
Nsj = length(Data);

nfig=0; 
mkdir figs 

%--------------------------------------------------------------------
% compare with surrogate data 
% 	  - frequence of choosing optimal path at each depth 
%--------------------------------------------------------------------
nfig=nfig+1; figure(nfig);clf;

Z = Data(1).Z;

Dn = [Data.dn];
Sn = [Data.sn];
An = [Data.an];

ax=mysubplot(3,12,.1,.1,.85,.85);
xxs = linspace(-.2,.2,Nsj);
for ss=1:6
	for dd=3:5
		axes(ax((dd-3)*12+ss*2-1)); axis off; 
		axes(ax((dd-3)*12+ss*2)); 
		j = 1:2^dd;
		[foo,ri] = sort(sum(Z.R(1:dd,j,dd-2,ss)));

		% compaute frequency of true sequence choice across all subjects 
		n = histc(An(Dn==dd & Sn==ss),j);
		freq(1:2^dd,ss,dd) = n(ri)/sum(n);

		% compaute frequency of surrogate sequence choice across all subjects 
		for mdl=1:nModls
			n = zeros(2^dd,1);
			for sj=1:Nsj
				ansurr = [SurrogateData(1).(models(mdl).name).an];
				ansurr = ansurr(Dn(:,sj)==dd & Sn(:,sj)==ss,:);
				n = n + histc(ansurr(:),j);
			end
			freqm(1:2^dd,ss,dd,mdl) = n(ri)/sum(n);
		end

		hon
		m1=nanmean(xa(j,:,ss,dd-2),2); 
		m2=nanmean(xb(j,:,ss,dd-2),2); 
		hof

		cla; 
		h=barh(m1);
		set(h,'facecolor',.7*ones(1,3));
		hon
		plot([m1-s1 m1+s1]',((1:2^dd)'*[1 1])','k')
		plot(m2,1:2^dd,'b','linewidth',2);
		hof
		clear tx;
		for k=1:2^dd
			tx{k} = num2str(Z.R(1:dd,ri(k),ddi,ss)');
		end
		set(gca,'ytick',j,'yticklabel',tx);
		set(gca,'fontsize',7);
		mytightaxes;
		ylim([0.5 2^dd+.5]);
		xl = xlim;
		xlim([0 xl(2)]);
		

		box on
		set(gca,'xticklabel',[]);
		if ~NLLO(dd-2,ss);set(gca,'xcolor','r','ycolor','r','linewidth',2);end
		if dd==3; title(['State ' num2str(ss)]);end
		if ss==1; ylabel(['Depth ' num2str(dd)]);end
	end
end

myfig(gcf,'figs/SurrogateDataPlots');

%--------------------------------------------------------------------
% compare with surrogate data 
% 	  - distribution over paths chosen 
%--------------------------------------------------------------------
nfig=nfig+1; figure(nfig);clf;

