function surrogateDataPlots(Data,models,SurrogateData,bestmodel)

nModls = length(models);
Nsj = length(Data);

nfig=0; 

mkdir figs 

%--------------------------------------------------------------------
% plot parameters of best model 
%--------------------------------------------------------------------

nfig=nfig+1; figure(nfig);clf;
for mdl=1:nModls
	R.(models(mdl).name) = load(sprintf('fitResults/%s',models(mdl).name));
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
myfig(gcf,'figs/Parameters');

%--------------------------------------------------------------------
% plot parameters of best model 
%--------------------------------------------------------------------

clear all; 
dt = '151223';

%whichfig = 1; % Figure 1 - parameter correlations and comparisons models 1-3 
%whichfig = 2; % Figure 2 - distributions of parameters for model 2 
whichfig = 3; % Figure 3 - stage 1 repetition probabilities 
%whichfig = 4; % Figure 4 - statistics on how many subjects fit well / poorly 
%whichfig = 5; % Figure 5 - Compare time course of generated choices to true data 
%whichfig = 6; % Figure 6 - Compare distributions of parameters, ask if their CIs don't overlap

modl=2; 

% list of models that can be fitted 
llfunc{1} = 'm2b2alr2';	 					% standard Daw et al. 2011 model 
llfunc{2} = 'bmfbmb2alr';					% mf/mb reparametrization 
llfunc{3} = '2bmfbmb2ar';					% mf/mb/lambda reparametrization 
llfunc{4} = 'bmfbmb2alcr';					% mf(0)/mf(1)/mb, two learning rates 
llfunc{5} = '2bmfbmb3ar';					% mf(0)/mf(1)/mb with three learning rates 
llfunc{6} = 'bmfbmb2alrb';					% mf/mb reparametrization - with bias 
llfunc{7} = 'bmfbmbalr';					% mf/mb reparametrization - one mf weight for both stages; one learning rate only
llfunc{8} = 'bmfalr';						% mf only 
llfunc{9} = 'bmbar';							% mb only at stage 1 

%NumParams = 8; 								% number of parmeters 
T = 201; 										% number of choices per subject 

Ti{1} = {'\beta_1','\beta_2','\alpha_1','\alpha_2','\lambda','\omega','\rho'};
Ti{2} = {'mb','mf','\beta_2','\alpha_1','\alpha_2','\lambda','\rho'};
Ti{3} = {'mb','mf0','mf1','\beta_2','\alpha_1','\alpha_2','\rho'};
Ti{4} = Ti{2}; 
Ti{5} = {'mb','mf0','mf1','\beta_2','\alpha_1','\alpha_2','\alpha_3','\rho'};
Ti{6} = {'mb','mf','\beta_2','\alpha_1','\alpha_2','\lambda','\rho','bias'};
Ti{7} = {'mb','mf','\alpha','\lambda','\rho'};
Ti{8} = {'mf','\alpha','\lambda','\rho'};
Ti{9} = {'mb','mf_2','\alpha','\rho'};

tTi{1} = {'log(\beta_1)','log(\beta_2)','\sigma^{-1}(\alpha_1)','\sigma^{-1}(\alpha_2)','\sigma^{-1}(\lambda)','\sigma^{-1}(\omega)','\rho'};
tTi{2} = {'log(mb)','log(mf)','log(\beta_2)','\sigma^{-1}(\alpha_1)','\sigma^{-1}(\alpha_2)','\sigma^{-1}(\lambda)','\rho'};
tTi{3} = {'log(mb)','log(mf0)','log(mf1)','log(\beta_2)','\sigma^{-1}(\alpha_1)','\sigma^{-1}(\alpha_2)','\rho'};
tTi{4} = tTi{2}; 
tTi{5} = {'log(mb)','log(mf0)','log(mf1)','log(\beta_2)','\sigma^{-1}(\alpha_1)','\sigma^{-1}(\alpha_2)','\sigma^{-1}(\alpha_3)','\rho'};
tTi{6} = {'log(mb)','log(mf)','log(\beta_2)','\sigma^{-1}(\alpha_1)','\sigma^{-1}(\alpha_2)','\sigma^{-1}(\lambda)','\rho','bias'};
tTi{7} = {'log(mb)','log(mf)','\sigma^{-1}(\alpha)','\sigma^{-1}(\lambda)','\rho'};
tTi{8} = {'log(mf)','\sigma^{-1}(\alpha)','\sigma^{-1}(\lambda)','\rho'};
tTi{9} = {'log(mb)','log(mf_2)','\sigma^{-1}(\alpha)','\rho'};

modelname{1} = 'Standard';
modelname{2} = 'mb/mf(\lambda)';
modelname{3} = 'mb/mf(0)/mf(1)';
modelname{4} = 'mb/mf(\lambda)_c';
modelname{5} = 'mb/mf(0)/mf(1)/3\alpha';
modelname{6} = 'mb/mf(0)/mf(1)+bias';
modelname{7} = 'mb/mf(\lambda)/1\alpha';
modelname{8} = 'mf(\lambda) only/1\alpha';
modelname{9} = 'mb only/1\alpha';

files ={'fit.m','emfit.m','addsimulations.m','runbatch'}; 
for k=1:length(llfunc); files{end+1} = ['ll' llfunc{k} '.m']; end
for k=1:length(llfunc); files{end+1} = ['gen' llfunc{k} '.m']; end


elseif whichfig==3 

	for proj=1:2
		loadstr=['mat/' dt '_TwostepDFG_P' num2str(proj) '_ll' llfunc{modl}];
		clear stats bf NumSubj Data Asurrog Ssurrog Rsurrog 
		eval(['load ' loadstr ' stats bf NumSubj NumParams Data Asurrog Ssurrog Rsurrog']);
		Nsj = NumSubj; 
		Nit = size(Asurrog,4);

		for sk=1:Nsj;
			a = Data(sk).A; 	
			i = find(~isnan(sum(a)));
			a = a(:,i);
			s = Data(sk).S(:,i); 	
			r = Data(sk).R(:,i); 	
			freqt = Data(sk).trans(i); 
			rep1 = a(1,2:end)==a(1,1:end-1); 
			rr(1,sk) = sum(rep1 & freqt(1:end-1)==1 & r(1:end-1)==1)/sum(freqt(1:end-1)==1 & r(1:end-1)==1);
			rr(2,sk) = sum(rep1 & freqt(1:end-1)==0 & r(1:end-1)==1)/sum(freqt(1:end-1)==0 & r(1:end-1)==1);
			rr(3,sk) = sum(rep1 & freqt(1:end-1)==1 & r(1:end-1)==0)/sum(freqt(1:end-1)==1 & r(1:end-1)==0);
			rr(4,sk) = sum(rep1 & freqt(1:end-1)==0 & r(1:end-1)==0)/sum(freqt(1:end-1)==0 & r(1:end-1)==0);
			for it=1:Nit
				a = Asurrog(:,:,sk,it); 
				i = ~isnan(sum(a));
				a = a(:,i);
				s = Ssurrog(:,:,sk,it); 
				r = Rsurrog(:,sk,it)';
				freqt = Data(sk).trans(i); 
				rep1 = a(1,2:end)==a(1,1:end-1); 
				rrs(1,sk,it) = sum(rep1 & freqt(1:end-1)==1 & r(1:end-1)==1)/sum(freqt(1:end-1)==1 & r(1:end-1)==1);
				rrs(2,sk,it) = sum(rep1 & freqt(1:end-1)==0 & r(1:end-1)==1)/sum(freqt(1:end-1)==0 & r(1:end-1)==1);
				rrs(3,sk,it) = sum(rep1 & freqt(1:end-1)==1 & r(1:end-1)==0)/sum(freqt(1:end-1)==1 & r(1:end-1)==0);
				rrs(4,sk,it) = sum(rep1 & freqt(1:end-1)==0 & r(1:end-1)==0)/sum(freqt(1:end-1)==0 & r(1:end-1)==0);
			end
		end

		mb = rr(1,:) - rr(2,:) - (rr(3,:) - rr(4,:));
		mf = rr(1,:) + rr(2,:) - (rr(3,:) + rr(4,:));
		mbs = rrs(1,:,:) - rrs(2,:,:) - (rrs(3,:,:) - rrs(4,:,:));
		mfs = rrs(1,:,:) + rrs(2,:,:) - (rrs(3,:,:) + rrs(4,:,:));
		clear m s 
		m(1,:) = mean(rr');
		m(2,:) = mean(mean(rrs,3),2);
		clear s
		s(1,:) = std(rr')/sqrt(Nsj);
		s(2,:) = std(mean(rrs,3),[],2)/sqrt(Nsj);

		xx = (1:4)'*ones(1,NumSubj) + .5*ones(4,1)*rand(1,NumSubj) -.25; 

		figure(3)

		subplot(2,7,(proj-1)*7+1);
		mybar(mean(rr'),.3);
		hon
		plot(xx,rr,'k.');
		hof
		xlim([.5 4.5]);
		title('True');
		set(gca,'xticklabel',{'Rew/Freq','Rew/Rare','Pun/Freq','Pun/Rare'},'xticklabelrotation',30);
		ylabel({['Project ' num2str(proj)],'Repeat probability'})

		subplot(2,7,(proj-1)*7+2);
		mybar(mean(mean(rrs,3),2),.7);
		hon
		plot(xx,mean(rrs,3),'k.');
		hof
		xlim([.5 4.5]);
		title('Surrogate');
		set(gca,'xticklabel',{'Rew/Freq','Rew/Rare','Pun/Freq','Pun/Rare'},'xticklabelrotation',30);

		subplot(2,7,(proj-1)*7+3);
		mybar(m',[.3;.7]);
		hon
		mydeb(0,m',s');
		hof
		xlim([.5 4.5]);
		title('Comparison');
		set(gca,'xticklabel',{'Rew/Freq','Rew/Rare','Pun/Freq','Pun/Rare'},'xticklabelrotation',30);

		subplot(2,7,(proj-1)*7+4);
		[c,p] = corr(mf',mean(mfs,3)'); 
		plot(mf',mean(mfs,3)','k.');
		title({'Reward effect',sprintf('c=%.2g p=%.2g',c,p)}); 
		xlabel('True data');
		ylabel('Surrogate data');
		mytightaxes; 

		subplot(2,7,(proj-1)*7+5);
		[c,p] = corr(mb',mean(mbs,3)'); 
		plot(mb',mean(mbs,3)','k.');
		title({'Reward x Frequency',sprintf('c=%.2g p=%.2g',c,p)}); 
		xlabel('True data');
		ylabel('Surrogate data');
		mytightaxes; 

		subplot(2,7,(proj-1)*7+6);
		[c,p] = corr(mf',stats.EMAP0','type','Spearman'); 
		mybar(c,.7);
		for k=1:NumParams; 
			if p(k)<.05; str='*'; 
			elseif p(k)<.05/7; str='**';
			else; str=''; 
			end
			h=text(k,c(k),sprintf('%s',str)); 
			set(h,'fontweight','bold','fontsize',18,'horizontalalignment','center');
		end
		set(gca,'xticklabel',tTi{modl},'xticklabelrotation',90);
		title({'Correlation w/','reward effect'}); 
		mytightaxes; 

		subplot(2,7,(proj-1)*7+7);
		[c,p] = corr(mb',stats.EMAP0','type','Spearman'); 
		mybar(c,.7);
		for k=1:NumParams; 
			if p(k)<.05; str='*'; 
			elseif p(k)<.05/7; str='**';
			else; str=''; 
			end
			h=text(k,c(k),sprintf('%s',str)); 
			set(h,'fontweight','bold','fontsize',18,'horizontalalignment','center');
		end
		set(gca,'xticklabel',tTi{modl},'xticklabelrotation',90);
		title({'Correlation w/','reward x freq int'}); 
		mytightaxes; 


	end

	myfig(gcf,['../figs/151224-SurrogateDataRepetitionEffects_' llfunc{modl}],'figcompare.m',files)

elseif whichfig==4

	figure(whichfig); clf; 

	for proj=1:2
		loadstr=['mat/' dt '_TwostepDFG_P' num2str(proj) '_ll' llfunc{modl}];
		clear stats bf NumSubj; eval(['load ' loadstr ' stats bf NumSubj Data Asurrog Ssurrog Rsurrog']);
		Nsj = NumSubj; 
		Nit = size(Asurrog,4);

		clear Nch 
		for k=1:Nsj; 
			Nch(k) = sum(sum(~isnan(Data(k).A)));
		end

		subplot(2,6,(proj-1)*3+(1:2))
			pc = exp(-stats.LL./Nch); 
			ps = 1-binocdf(pc.*Nch,Nch,.5);
			i = 1:NumSubj ; 
			j = ps<.05; 
			CapturedBetterThanChance = j; 
			plot(i,pc,'o');
			hon
			plot(i(~j),pc(~j),'r*');
			plot([1 NumSubj],median(pc)*[1 1],'k--','linewidth',2);
			hof
			axis tight 
			ylim([.45 .95])
			ylabel({'Average choice probability',['afforded by model' modelname{modl}]})
			title({['P' num2str(proj)],['#> .05: ' num2str(sum(j)) ' (' sprintf('%.2g',sum(j)/Nsj*100) '%)']})

		subplot(2,12,(proj-1)*6+5)
			xx = linspace(.4,1,40); 
			nn = histc(pc,xx);
			nn2 = histc(pc(~j),xx);
			barh(xx,nn/sum(nn));
			hon
			h=barh(xx,nn2/sum(nn));
			set(h,'facecolor','r'); 
			hof
			axis tight 
			set(gca,'xticklabel',[],'yticklabel',[]);
			ylim([.4 1])

		subplot(2,6,(proj-1)*3+(1:2)+6)
			plot(Nch,'o');
			hon
			plot(i(~j),Nch(~j),'r.')
			hof
			ylabel('Number of valid choices')
			xlabel('Subject')
			axis tight 
			title(['Valid choices (' sprintf('%.3g',mean(Nch)/402*100) '%)'])

	end

	myfig(gcf,['../figs/151224-IndividualFits_' llfunc{modl}],'figcompare.m',files)


elseif whichfig==5

		figure(5);clf; 
		col = jet(length(llfunc));

		for proj=1:2

			clear m n s ms ns ss; 
			modli = [1:3 6:9]; 1:length(llfunc); 

			for modl=modli; 
				loadstr=['mat/' dt '_TwostepDFG_P' num2str(proj) '_ll' llfunc{modl}]
				clear stats bf NumSubj; eval(['load ' loadstr ' stats NumSubj Data Asurrog Ssurrog Rsurrog']);
				Nsj = NumSubj; 
				Nit = size(Asurrog,4);

				AA = NaN*ones(201,Nsj,2);
				SS = NaN*ones(201,Nsj,2);
				for sk=1:Nsj;
					a = Data(sk).A; 	
					s = Data(sk).S; 	
					Nch = sum(~isnan(a(:)));
					pch(sk) = exp(-stats.LL(sk)/Nch);
					%sjfit = 1-binocdf(Nch*pch,Nch,.5)<.05;
					i = 1:length(a);
					AA(i,sk,:) = a'; 
					SS(i,sk,:) = s(:,i)'; 
				end
				mLL(modl,proj,1) = mean(stats.LL); 
				mLL(modl,proj,2) = mean(pch); 

				if modl==modli(1)

					n=sum(~isnan(AA(:,:,1)),2);
					m=sq(nansum( AA(:,:,1)==1    ,2));
					s=sq(nansum((AA(:,:,1)==1).^2,2));

					for k=1:2
						n(:,k+1)=    nansum(SS(:,:,2)==k+1,2);
						m(:,k+1)=sq(nansum( AA(:,:,2)==1&SS(:,:,2)==k+1    ,2));
						s(:,k+1)=sq(nansum((AA(:,:,2)==1&SS(:,:,2)==k+1).^2,2));
					end

					m=m./n; 
					s=s./n; 
					s=1.96*sqrt(s-m.^2)./sqrt(n);
					subplot(3,2,proj)
						h=fill([1:201 201:-1:1]',[m(:,1)+s(:,1);m(end:-1:1,1)-s(end:-1:1,1)],[.8 .8 1],'edgecolor',[.8 .8 .8]); 
						hon; 
						plot(m(:,1),'k.');
						mytightaxes
						ylabel(['P(action 1 in State ' num2str(1) ')']);
						%ylim([.25 .9])
						title(['P' num2str(proj)]);
						for k=1:2
						subplot(3,2,proj+2+(k-1)*2)
							h=fill([1:201 201:-1:1]',[m(:,k+1)+s(:,k+1);m(end:-1:1,k+1)-s(end:-1:1,k+1)],[.8 .8 1],'edgecolor',[.8 .8 .8]); 
							hon; 
							plot(m(:,k+1),'k.');
							mytightaxes
							ylabel(['P(action 1 in State ' num2str(k+1) ')']);
							if k==2; xlabel('Trial');end
							%ylim([.25 .9])

						end
				end

				ms(:,1)=nansum(nansum( Asurrog(1,:,:,:)==1    ,4),3)';
				ss(:,1)=nansum(nansum((Asurrog(1,:,:,:)==1).^2,4),3)';
				if any(isnan(Asurrog(:)));error;end
				ns(:,1)=NumSubj*Nit*ones(T,1);
				for k=1:2
					ms(:,k+1)=nansum(nansum( Asurrog(2,:,:,:)==1&Ssurrog(2,:,:,:)==k    ,4),3);
					ss(:,k+1)=nansum(nansum((Asurrog(2,:,:,:)==1&Ssurrog(2,:,:,:)==k).^2,4),3);
					ns(:,k+1)=nansum(nansum(                     Ssurrog(2,:,:,:)==k     ,4),3);
				end
				ms = ms./ns; 
				ss = 1.96*(ss./ns-ms.^2);%./sqrt(ns); 


				subplot(3,2,proj)
				hh(modl) = plot(ms(:,1),'linewidth',1,'color',col(modl,:));
				%plot(ms(:,1)+ss(:,1),'--','linewidth',1,'color',col(modl,:));
				%plot(ms(:,1)-ss(:,1),'--','linewidth',1,'color',col(modl,:));
				if modl==length(llfunc); hof; end
				for k=1:2
					subplot(3,2,proj+2+(k-1)*2)
						plot(ms(:,k+1),'linewidth',1,'color',col(modl,:));
				%		plot(ms(:,k+1)+ss(:,k+1),'--','linewidth',1,'color',col(modl,:));
				%		plot(ms(:,k+1)-ss(:,k+1),'--','linewidth',1,'color',col(modl,:));
						if modl==length(llfunc); hof; end
				end

			end
		end
		%legend(hh([1 2 6]),llfunc{[1,2,6]})
		legend(hh(modli),modelname(modli))
	myfig(gcf,'../figs/151223-LearningCurves','figcompare.m',files); 


elseif whichfig==6 

	paramtransform = ['pt' llfunc{modl} ];
	paramtransform = str2func(paramtransform);
	p = [2 3]; 
	p = [4 5]; 
	p = [1 2]; 
		
	for proj=1:2
		subplot(1,2,proj)
		loadstr=['mat/' dt '_TwostepDFG_P' num2str(proj) '_ll' llfunc{modl}];
		eval(['load ' loadstr ' stats']);
		E = stats.EMAP0; 
		V = 1.96*sqrt(stats.VMAP0); 

		plot(ones(2,1)*E(p(1),:),[E(p(2),:)-V(p(2),:); E(p(2),:)+V(p(2),:)],'color',[.7 .7 1]);
		hon
		plot([E(p(1),:)-V(p(1),:); E(p(1),:)+V(p(1),:)],ones(2,1)*E(p(2),:),'color',[.7 .7 1]);
		plot(E(p(1),:),E(p(2),:),'k.','markersize',15); 


		i = 	(E(p(2),:) < E(p(1),:)-V(p(1),:)) | (E(p(2),:) > E(p(1),:)+V(p(1),:)) & ...
		  		(E(p(1),:) < E(p(2),:)-V(p(2),:)) | (E(p(1),:) > E(p(2),:)+V(p(2),:)); 
		i = 	(E(p(2),:)+V(p(2),:) < E(p(1),:)-V(p(1),:)) | (E(p(2),:)-V(p(2),:) > E(p(1),:)+V(p(1),:)); 

		plot(E(p(1),i),E(p(2),i),'r.','markersize',20); 
		plot(ones(2,1)*E(p(1),i),[E(p(2),i)-V(p(2),i); E(p(2),i)+V(p(2),i)],'r');
		hon
		plot([E(p(1),i)-V(p(1),i); E(p(1),i)+V(p(1),i)],ones(2,1)*E(p(2),i),'r');
		plot([-6 6],[-6 6],'r','linewidth',2);
		hof

		sum(i)/length(i)
		mytightaxes; 
		xlabel(tTi{modl}(p(1)))
		ylabel(tTi{modl}(p(2)))
		title(sprintf('P%i: %.2g%% not different',proj,100-sum(i)/length(i)*100))
	end
	myfig(gcf,'../figs/151223-LearningRateDifferences','figcompare.m',files)
	%myfig(gcf,'../figs/151223-MFBetaDifferences','figcompare.m',files)

end

