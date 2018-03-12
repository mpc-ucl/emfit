function surrogateDataPlots(Data,models,SurrogateData,bestmodel)

nModls = length(models);
Nsj = length(Data);

nfig=0; 

mkdir figs 

%--------------------------------------------------------------------
% compare with surrogate data 
%--------------------------------------------------------------------
nfig=nfig+1; figure(nfig);clf;

T = length(Data(1).a); 
ns = zeros(size([Data.a],1)/2,2);
nns = zeros(size([Data.a],1)/2,2,nModls);

% determine block length for blocked response bias plots 
blocklength = floor(length(Data(1).a)/3); 
if blocklength*3 ~=length(Data(1).a); 
	fprintf('warning: block length %i but data length %i - throwing away last few trials!\n',blocklength,length(Data(1).a));
end

wdw = 20; 

for sj=1:Nsj
	a = Data(sj).a; 
	s = Data(sj).s; 

	% switch so action 1 is correct response to stim 1 
	if Data(sj).I(1,2)==1; 
		a = 3-a; 
	end

	% switch so stim 1 is more rewarded 
	prc = Data(sj).prc; 
	switchstim=0; 
	if prc(2)>prc(1); 
		s = 3-s; 
		a = 3-a; 
		prc = prc(2:-1:1);
		switchstim=1; 
	end

	for ss=1:2
		i = s==ss; 
		as(1:sum(i),ss,sj) = a(i)==1;
		ns(1:sum(i),ss) = ns(1:sum(i),ss)+1;
	end

	% response bias in blocks 
	for bl=1:3
		i = (1:blocklength)+(bl-1)*blocklength; 
		bias(bl,sj) = 0.5 * log(...
			(nansum(a(i)==1 & s(i)==1)+.5) *...
			(nansum(a(i)==1 & s(i)==2)+.5) /...
			(nansum(a(i)==2 & s(i)==1)+.5) /...
			(nansum(a(i)==2 & s(i)==2)+.5) );
	end

	% response bias - sliding window 
	j=0;
	for bl=1:wdw/2:(length(Data(1).a)-wdw)
		j=j+1;
		i = (1:wdw)+(bl-1); 
		wdwi(j) = i(1);
		bias_sl(j,sj) = 0.5 * log(...
			(nansum(a(i)==1 & s(i)==1)+.5) *...
			(nansum(a(i)==1 & s(i)==2)+.5) /...
			(nansum(a(i)==2 & s(i)==1)+.5) /...
			(nansum(a(i)==2 & s(i)==2)+.5) );
	end


	for mdl=1:nModls;
		a = [SurrogateData(sj).(models(mdl).name).a];

		% switch so action 1 is correct response to stim 1 
		if Data(sj).I(1,2)==1; 
			a = 3-a; 
		end
		% switched stimuli so action 1 is more rewarded? 
		if switchstim
			a = 3-a; 
		end

		b = mean(a==1,2);
		for ss=1:2
			i = s==ss; 
			bs(1:sum(i),ss,mdl,sj) = b(i);
			nns(1:sum(i),ss,mdl) = nns(1:sum(i),ss,mdl)+1;
		end

		% response bias in blocks 
		for bl=1:3
			i = (1:blocklength)+(bl-1)*blocklength; 
			for nsample=1:size(a,2);
				biass(bl,mdl,sj,nsample) = 0.5 * log(...
					(nansum(a(i,nsample)==1 & s(i)==1)+.5) *...
					(nansum(a(i,nsample)==1 & s(i)==2)+.5) /...
					(nansum(a(i,nsample)==2 & s(i)==1)+.5) /...
					(nansum(a(i,nsample)==2 & s(i)==2)+.5) );
			end
		end
		% response bias - sliding window 
		j=0;
		for bl=1:wdw/2:(length(Data(1).a)-wdw)
			j=j+1;
			i = (1:wdw)+(bl-1); 
			for nsample=1:size(a,2)
				biass_sl(j,mdl,sj,nsample) = 0.5 * log(...
					(nansum(a(i,nsample)==1 & s(i)==1)+.5) *...
					(nansum(a(i,nsample)==1 & s(i)==2)+.5) /...
					(nansum(a(i,nsample)==2 & s(i)==1)+.5) /...
					(nansum(a(i,nsample)==2 & s(i)==2)+.5) );
			end
		end
	end

end
mas = sum(as,3)./ns;
mbs = sum(bs,4)./nns;
mbias = mean(bias,2);
sbias = std(bias,[],2)/sqrt(Nsj);
mbiass = mean(mean(biass,4),3);
sbiass = std(mean(biass,4),[],3);

mbias_sl = mean(bias_sl,2);
mbiass_sl = mean(mean(biass_sl,4),3);

Ti = {'Rich stimulus','Lean stimulus'}; 

subplot(1,3,1)
	mybar(mbias,.7);
	hon
	xx = [1:3]'*ones(1,nModls) + ones(3,1)*linspace(-.3,.3,nModls);
	plot(xx,mbiass,'.-','markersize',15,'linewidth',1)
	hof
	xlim([.5 3.5]);
	ylabel('Response bias');
	set(gca,'xtick',1:3);
	xlabel(sprintf('Block of %i trials',blocklength));
		
for ss=1:2
	subplot(2,3,1+ss);
		plot(mas(:,ss),'k','linewidth',3);
		hon
		plot(sq(mbs(:,ss,:)),'linewidth',2)
		hof
		ylim([0 1]);
		title(Ti{ss});
		xlabel('Trial');
		if ss==1; ylabel('P(rich rsponse)');end
end

subplot(2,3,5:6);
	plot(mbias_sl,'k','linewidth',3);
	hon
	plot(mbiass_sl,'linewidth',2)
	hof
	set(gca,'xtick',1:2:length(wdwi),'xticklabel',wdwi(1:2:end))
	title('Response bias over sliding window');
	xlabel('Trial');

le = {models.name}; 
le = {'Data',le{:}};
legend(le,'location','best'); 

myfig(gcf,'figs/SurrogateDataPlots');


