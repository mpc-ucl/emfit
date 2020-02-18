function surrogateDataPlots(Data,models,SurrogateData,bestmodel,fitResults); 

nModls = length(models);
Nsj = length(Data);

nfig=0; 
mkdir figs 

%--------------------------------------------------------------------
% compare generated reward-dependent choices to real ones 
% do so across all subjects, and separately for those fit best by flat and non-flat models 
%--------------------------------------------------------------------

for sj=1:Nsj
	a = Data(sj).a';
	for rew = 3:7
		i = rew==Data(sj).rew;
		probHE(rew-2,sj)	= sum(a(i,:)==2)/sum(i);
	end
end

% get choice probability as function of reward in surrogate data, averaging over
% the samples
for mdl=1:nModls
	for sj=1:Nsj
		Asurr = [SurrogateData(sj).(models(mdl).name).a]';
		for rew = 3:7
			i = rew==Data(sj).rew;
			probHESurr(rew-2,sj,mdl)	= mean(sum(Asurr(i,:)==2)/sum(i));
		end
	end
end

rews = [3:7];

nfig=nfig+1; clf;
subplot(121)
	m = mean(probHE,2);
	ms = squeeze(mean(probHESurr,2));
	plot(rews,m,'k.-','linewidth',2,'markersize',20);
	hold on
	plot(rews,ms);
	hold off
	xlabel({'Reward for','high effort option'});
	ylabel('% chosen high effort');
	set(gca,'fontsize',14);
	mytightaxes; 
    
    
%--------------------------------------------------------------------
% compare generated decision times to real ones 
%--------------------------------------------------------------------

for sj=1:Nsj
	dt = Data(sj).decisiontime;
	for rew = 3:7
		i = rew==Data(sj).rew;
		decisionTimes(rew-2,sj)	= sum(dt(i))/sum(i);
	end
end

% get choice probability as function of reward in surrogate data, averaging over
% the samples
for mdl=1:nModls
   nSamples = size(SurrogateData(sj).(models(mdl).name),2);
	for sj=1:Nsj
		for k = 1:nSamples
            AsurrTime = SurrogateData(sj).(models(mdl).name)(k).simTime(:);
            validTrials = AsurrTime>0.7;
            for rew = 3:7
                i = rew==Data(sj).rew;
					 i = i&validTrials; 
                decisionTimeSurrPerGen(rew-2,sj,mdl,k)	= mean(AsurrTime(i));
            end
        end
        decisionTimesSurr(:,sj,mdl)	= mean(decisionTimeSurrPerGen(:,sj,mdl,:),4);
    end
end

rews = [3:7];

subplot(122)
	m = nanmean(decisionTimes,2);
	ms = squeeze(nanmean(decisionTimesSurr,2));
	plot(rews,m,'k.-','linewidth',2,'markersize',20);
	hold on
	plot(rews,ms);
	hold off
	xlabel({'Reward for','high effort option'});
	ylabel('decision time (sec)');
	set(gca,'fontsize',14);
	h = legend({'Data',models.name},'location','best');
	set(h,'fontsize',10);
	mytightaxes; 
    
    
%% save figure

myfig(gcf,sprintf('%s/figs/SurrogateChoiceDecisionTime',fitResults));
