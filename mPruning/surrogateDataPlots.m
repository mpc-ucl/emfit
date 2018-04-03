function surrogateDataPlots(Data,models,SurrogateData,bestmodel,resultsDir)

nModls = length(models);
Nsj = length(Data);

nfig=0; 
mkdir figs 

%--------------------------------------------------------------------
% compare with surrogate data 
% 	  - distribution over paths chosen 
%--------------------------------------------------------------------
nfig=nfig+1; figure(nfig);clf;

Z = Data(1).Z;

Dn = [Data.dn];
Sn = [Data.sn];
An = [Data.an];

ax=mysubplot(4,12,.1,.1,.85,.85);
xxs = linspace(-.2,.2,Nsj);
for ss=1:6
	for dd=3:6
		axes(ax((dd-3)*12+ss*2-1)); axis off; 
		axes(ax((dd-3)*12+ss*2)); 
		j = 1:2^dd;
		[foo,ri] = sort(sum(Z.R(1:dd,j,dd-2,ss)));

		% compute frequency of true sequence choice across all subjects 
		n = histc(An(Dn==dd & Sn==ss),j);
		freq(1:2^dd,ss,dd) = n(ri)/sum(n);

		% compute frequency of surrogate sequence choice across all subjects 
		for mdl=1:nModls
			n = zeros(2^dd,1);
			for sj=1:Nsj
				ansurr = [SurrogateData(sj).(models(mdl).name).an];
				ansurr = ansurr(Dn(:,sj)==dd & Sn(:,sj)==ss,:);
				n = n + histc(ansurr(:),j);
			end
			freqm(1:2^dd,ss,dd,mdl) = n(ri)/sum(n);
		end

		cla; 
		h=barh(freq(1:2^dd,ss,dd));
		set(h,'facecolor',.7*ones(1,3));
		hon
		plot(sq(freqm(1:2^dd,ss,dd,:)),1:2^dd,'linewidth',2);
		hof

		% add labels with rewards along each path 
		clear tx;
		for k=1:2^dd
			tx{k} = num2str(Z.R(1:dd,ri(k),dd-2,ss)');
		end
		set(gca,'ytick',j,'yticklabel',tx);
		set(gca,'fontsize',7);
		mytightaxes;
		ylim([0.5 2^dd+.5]);
		xl = xlim;
		xlim([0 xl(2)]);

		% colour state/depth combinations where optimum involes large loss
		% transition red 
		box on
		set(gca,'xticklabel',[]);
		best = ri(end);
		rbest = Z.R(1:dd,best,dd-2,ss);
		if any(rbest==Z.r3);
			set(gca,'xcolor','r','ycolor','r','linewidth',2);
		end
		if dd==3; title(['State ' num2str(ss)]);end
		if ss==1; ylabel(['Depth ' num2str(dd)]);end
	end
end

myfig(gcf,'figs/SurrogateDataPlots');

%--------------------------------------------------------------------
% compare with surrogate data 
% 	  - frequence of choosing optimal path at each depth 
%--------------------------------------------------------------------

clear n  best; 
for dd=3:6
	j = 1:2^dd;
	for ss=1:6
		rews = Z.R(1:dd,j,dd-2,ss);
		[foo,ri] = sort(sum(rews));
		best(ss,dd) = ri(end);
		if any(rews(:,ri(end))==Z.r3); 
			LLO(ss,dd)=1;
		end
	end
end

for sj=1:Nsj 
	an = Data(sj).an;
	sn = Data(sj).sn;
	dn = Data(sj).dn;
	for dd=3:6
		for ss=1:6
			ncorr(ss,dd) = sum(an(sn==ss & dn==dd)==best(ss,dd));
			n(ss,dd) = sum(sn==ss & dn==dd); 
		end
	end
	for dd=3:6
		pcorr(dd-2,1,sj) = sum(ncorr(LLO(:,dd)==0,dd)) / sum(n(LLO(:,dd)==0,dd));
		pcorr(dd-2,2,sj) = sum(ncorr(LLO(:,dd)==1,dd)) / sum(n(LLO(:,dd)==1,dd));
	end

	for mdl=1:nModls
		ansurr = [SurrogateData(sj).(models(mdl).name).an];
		for dd=3:6
			for ss=1:6
				ncorrs(ss,dd) = sum(sum(ansurr(sn==ss & dn==dd,:)==best(ss,dd)));
			end
		end
		for dd=3:6
			pcorrs(dd-2,1,mdl,sj) = sum(ncorrs(LLO(:,dd)==0,dd)) / sum(n(LLO(:,dd)==0,dd)) / size(ansurr,2);
			pcorrs(dd-2,2,mdl,sj) = sum(ncorrs(LLO(:,dd)==1,dd)) / sum(n(LLO(:,dd)==1,dd)) / size(ansurr,2);
		end
	end
end

m = nanmean(pcorr,3);
s = nanstd(pcorr,[],3)/sqrt(Nsj);
ms = nanmean(pcorrs,4);

nfig=nfig+1; figure(nfig);clf;
	h1=bar(m);
	hon
	mydeb(0,m,s);
	h2=plot(sq(ms(:,1,:)));
	plot(sq(ms(:,2,:)));
	hof
	legend(h1,'No large loss','Large loss');
	xlabel('Fraction optimal choice');
	ylabel('Depth');
	set(gca,'xticklabel',[3:6]);


