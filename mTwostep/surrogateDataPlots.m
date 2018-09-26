function surrogateDataPlots(Data,models,SurrogateData,bestmodel,resultsDir)

nModls = length(models);
Nsj = length(Data);

nfig=get(gcf,'Number');
mkdir figs 


%--------------------------------------------------------------------
% stage 1 repetition probabilities 
%--------------------------------------------------------------------

nfig=nfig+1; figure(nfig);clf;

nSample=length(SurrogateData(1).(models(1).name));

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
	for mdl=1:nModls
		for it=1:nSample
			a = SurrogateData(sk).(models(mdl).name)(it).A;
			i = ~isnan(sum(a));
			a = a(:,i);
			s = SurrogateData(sk).(models(mdl).name)(it).S(i);
			r = SurrogateData(sk).(models(mdl).name)(it).R(i);
			freqt = Data(sk).trans(i); 
			rep1 = a(1,2:end)==a(1,1:end-1); 
			rrs(1,mdl,sk,it) = sum(rep1 & freqt(1:end-1)==1 & r(1:end-1)==1)/sum(freqt(1:end-1)==1 & r(1:end-1)==1);
			rrs(2,mdl,sk,it) = sum(rep1 & freqt(1:end-1)==0 & r(1:end-1)==1)/sum(freqt(1:end-1)==0 & r(1:end-1)==1);
			rrs(3,mdl,sk,it) = sum(rep1 & freqt(1:end-1)==1 & r(1:end-1)==0)/sum(freqt(1:end-1)==1 & r(1:end-1)==0);
			rrs(4,mdl,sk,it) = sum(rep1 & freqt(1:end-1)==0 & r(1:end-1)==0)/sum(freqt(1:end-1)==0 & r(1:end-1)==0);
		end
	end
end

mb = rr(1,:) - rr(2,:) - (rr(3,:) - rr(4,:));
mf = rr(1,:) + rr(2,:) - (rr(3,:) + rr(4,:));
mbs = rrs(1,:,:,:) - rrs(2,:,:,:) - (rrs(3,:,:,:) - rrs(4,:,:,:));
mfs = rrs(1,:,:,:) + rrs(2,:,:,:) - (rrs(3,:,:,:) + rrs(4,:,:,:));

clear m s 
m(1,:) = nanmean(rr');
m(1+(1:nModls),:) = nanmean(nanmean(rrs,4),3)';
clear s
s(1,:) = nanstd(rr')/sqrt(Nsj);
s(1+(1:nModls),:) = nanstd(nanmean(rrs,4),[],3)'/sqrt(Nsj);

xx = (1:4)'*ones(1,Nsj) + .5*ones(4,1)*rand(1,Nsj) -.25; 


subplot(2,4,1:4);
h=bar(m'); 
hon
mydeb(0,m',s');
hof
l{1} = 'Data';
l(1+(1:nModls)) = {models.name};
legend(l);
xlim([.5 4.5]);
ylim([0 1]);
title('Comparison');
ylabel('Level 1 repeat probability');
set(gca,'xticklabel',{'Rew/Freq','Rew/Rare','Pun/Freq','Pun/Rare'},'xticklabelrotation',30);


subplot(2,4,5);
[c,p] = corr(mf',sq(nanmean(mfs(:,bestmodel,:,:),4))); 
plot(mf',sq(nanmean(mfs(:,bestmodel,:,:),4)),'k.');
title({'Reward effect',sprintf('c=%.2g p=%.2g',c,p)}); 
xlabel('True data');
ylabel('Surrogate data');
mytightaxes; 

subplot(2,4,6);
[c,p] = corr(mb',sq(mean(mbs(:,bestmodel,:,:),4))); 
plot(mb',sq(mean(mbs(:,bestmodel,:,:),4)),'k.');
title({'Reward x Frequency',sprintf('c=%.2g p=%.2g',c,p)}); 
xlabel('True data');
ylabel('Surrogate data');
mytightaxes; 

subplot(2,4,7);

for mdl=1:nModls
	R.(models(mdl).name) = load(sprintf('%s/%s',resultsDir,models(mdl).name));
end

[c,p] = corr(mf',R.(models(bestmodel).name).stats.EMAP0','type','Spearman'); 
mybar(c,.7);
for k=1:models(bestmodel).npar
	if p(k)<.05; str='*'; 
	elseif p(k)<.05/7; str='**';
	else; str=''; 
	end
	h=text(k,c(k),sprintf('%s',str)); 
	set(h,'fontweight','bold','fontsize',18,'horizontalalignment','center');
end
set(gca,'xticklabel',models(bestmodel).parnames,'xticklabelrotation',90);
title({'Correlation w/','reward effect'}); 
mytightaxes; 

subplot(2,4,8);
[c,p] = corr(mb',R.(models(bestmodel).name).stats.EMAP0','type','Spearman'); 
mybar(c,.7);
for k=1:models(bestmodel).npar
	if p(k)<.05; str='*'; 
	elseif p(k)<.05/7; str='**';
	else; str=''; 
	end
	h=text(k,c(k),sprintf('%s',str)); 
	set(h,'fontweight','bold','fontsize',18,'horizontalalignment','center');
end
set(gca,'xticklabel',models(bestmodel).parnames,'xticklabelrotation',90);
title({'Correlation w/','reward x freq int'}); 
mytightaxes; 

myfig(gcf,[resultsDir filesep 'figs/TwostepSurrogateDataRepetitionEffects']);

% elseif whichfig==5
% %whichfig = 5; % Figure 5 - Compare time course of generated choices to true data 
% 
% 		figure(5);clf; 
% 		col = jet(length(llfunc));
% 
% 		for proj=1:2
% 
% 			clear m n s ms ns ss; 
% 			modli = [1:3 6:9]; 1:length(llfunc); 
% 
% 			for modl=modli; 
% 				loadstr=['mat/' dt '_TwostepDFG_P' num2str(proj) '_ll' llfunc{modl}]
% 				clear stats bf Nsj; eval(['load ' loadstr ' stats Nsj Data Asurrog Ssurrog Rsurrog']);
% 				Nsj = Nsj; 
% 				Nit = size(Asurrog,4);
% 
% 				AA = NaN*ones(201,Nsj,2);
% 				SS = NaN*ones(201,Nsj,2);
% 				for sk=1:Nsj;
% 					a = Data(sk).A; 	
% 					s = Data(sk).S; 	
% 					Nch = sum(~isnan(a(:)));
% 					pch(sk) = exp(-stats.LL(sk)/Nch);
% 					%sjfit = 1-binocdf(Nch*pch,Nch,.5)<.05;
% 					i = 1:length(a);
% 					AA(i,sk,:) = a'; 
% 					SS(i,sk,:) = s(:,i)'; 
% 				end
% 				mLL(modl,proj,1) = mean(stats.LL); 
% 				mLL(modl,proj,2) = mean(pch); 
% 
% 				if modl==modli(1)
% 
% 					n=sum(~isnan(AA(:,:,1)),2);
% 					m=sq(nansum( AA(:,:,1)==1    ,2));
% 					s=sq(nansum((AA(:,:,1)==1).^2,2));
% 
% 					for k=1:2
% 						n(:,k+1)=    nansum(SS(:,:,2)==k+1,2);
% 						m(:,k+1)=sq(nansum( AA(:,:,2)==1&SS(:,:,2)==k+1    ,2));
% 						s(:,k+1)=sq(nansum((AA(:,:,2)==1&SS(:,:,2)==k+1).^2,2));
% 					end
% 
% 					m=m./n; 
% 					s=s./n; 
% 					s=1.96*sqrt(s-m.^2)./sqrt(n);
% 					subplot(3,2,proj)
% 						h=fill([1:201 201:-1:1]',[m(:,1)+s(:,1);m(end:-1:1,1)-s(end:-1:1,1)],[.8 .8 1],'edgecolor',[.8 .8 .8]); 
% 						hon; 
% 						plot(m(:,1),'k.');
% 						mytightaxes
% 						ylabel(['P(action 1 in State ' num2str(1) ')']);
% 						%ylim([.25 .9])
% 						title(['P' num2str(proj)]);
% 						for k=1:2
% 						subplot(3,2,proj+2+(k-1)*2)
% 							h=fill([1:201 201:-1:1]',[m(:,k+1)+s(:,k+1);m(end:-1:1,k+1)-s(end:-1:1,k+1)],[.8 .8 1],'edgecolor',[.8 .8 .8]); 
% 							hon; 
% 							plot(m(:,k+1),'k.');
% 							mytightaxes
% 							ylabel(['P(action 1 in State ' num2str(k+1) ')']);
% 							if k==2; xlabel('Trial');end
% 							%ylim([.25 .9])
% 
% 						end
% 				end
% 
% 				ms(:,1)=nansum(nansum( Asurrog(1,:,:,:)==1    ,4),3)';
% 				ss(:,1)=nansum(nansum((Asurrog(1,:,:,:)==1).^2,4),3)';
% 				if any(isnan(Asurrog(:)));error;end
% 				ns(:,1)=Nsj*Nit*ones(T,1);
% 				for k=1:2
% 					ms(:,k+1)=nansum(nansum( Asurrog(2,:,:,:)==1&Ssurrog(2,:,:,:)==k    ,4),3);
% 					ss(:,k+1)=nansum(nansum((Asurrog(2,:,:,:)==1&Ssurrog(2,:,:,:)==k).^2,4),3);
% 					ns(:,k+1)=nansum(nansum(                     Ssurrog(2,:,:,:)==k     ,4),3);
% 				end
% 				ms = ms./ns; 
% 				ss = 1.96*(ss./ns-ms.^2);%./sqrt(ns); 
% 
% 
% 				subplot(3,2,proj)
% 				hh(modl) = plot(ms(:,1),'linewidth',1,'color',col(modl,:));
% 				%plot(ms(:,1)+ss(:,1),'--','linewidth',1,'color',col(modl,:));
% 				%plot(ms(:,1)-ss(:,1),'--','linewidth',1,'color',col(modl,:));
% 				if modl==length(llfunc); hof; end
% 				for k=1:2
% 					subplot(3,2,proj+2+(k-1)*2)
% 						plot(ms(:,k+1),'linewidth',1,'color',col(modl,:));
% 				%		plot(ms(:,k+1)+ss(:,k+1),'--','linewidth',1,'color',col(modl,:));
% 				%		plot(ms(:,k+1)-ss(:,k+1),'--','linewidth',1,'color',col(modl,:));
% 						if modl==length(llfunc); hof; end
% 				end
% 
% 			end
% 		end
% 		%legend(hh([1 2 6]),llfunc{[1,2,6]})
% 		legend(hh(modli),modelname(modli))
% 	myfig(gcf,'../figs/151223-LearningCurves','figcompare.m',files); 
% 
% 
% elseif whichfig==6 
% 
% 	paramtransform = ['pt' llfunc{modl} ];
% 	paramtransform = str2func(paramtransform);
% 	p = [2 3]; 
% 	p = [4 5]; 
% 	p = [1 2]; 
% 		
% 	for proj=1:2
% 		subplot(1,2,proj)
% 		loadstr=['mat/' dt '_TwostepDFG_P' num2str(proj) '_ll' llfunc{modl}];
% 		eval(['load ' loadstr ' stats']);
% 		E = stats.EMAP0; 
% 		V = 1.96*sqrt(stats.VMAP0); 
% 
% 		plot(ones(2,1)*E(p(1),:),[E(p(2),:)-V(p(2),:); E(p(2),:)+V(p(2),:)],'color',[.7 .7 1]);
% 		hon
% 		plot([E(p(1),:)-V(p(1),:); E(p(1),:)+V(p(1),:)],ones(2,1)*E(p(2),:),'color',[.7 .7 1]);
% 		plot(E(p(1),:),E(p(2),:),'k.','markersize',15); 
% 
% 
% 		i = 	(E(p(2),:) < E(p(1),:)-V(p(1),:)) | (E(p(2),:) > E(p(1),:)+V(p(1),:)) & ...
% 		  		(E(p(1),:) < E(p(2),:)-V(p(2),:)) | (E(p(1),:) > E(p(2),:)+V(p(2),:)); 
% 		i = 	(E(p(2),:)+V(p(2),:) < E(p(1),:)-V(p(1),:)) | (E(p(2),:)-V(p(2),:) > E(p(1),:)+V(p(1),:)); 
% 
% 		plot(E(p(1),i),E(p(2),i),'r.','markersize',20); 
% 		plot(ones(2,1)*E(p(1),i),[E(p(2),i)-V(p(2),i); E(p(2),i)+V(p(2),i)],'r');
% 		hon
% 		plot([E(p(1),i)-V(p(1),i); E(p(1),i)+V(p(1),i)],ones(2,1)*E(p(2),i),'r');
% 		plot([-6 6],[-6 6],'r','linewidth',2);
% 		hof
% 
% 		sum(i)/length(i)
% 		mytightaxes; 
% 		xlabel(tTi{modl}(p(1)))
% 		ylabel(tTi{modl}(p(2)))
% 		title(sprintf('P%i: %.2g%% not different',proj,100-sum(i)/length(i)*100))
% 	end
% 	myfig(gcf,'../figs/151223-LearningRateDifferences','figcompare.m',files)
% 	%myfig(gcf,'../figs/151223-MFBetaDifferences','figcompare.m',files)
% 
% end
% 
