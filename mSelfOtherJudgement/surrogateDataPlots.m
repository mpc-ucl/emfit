function surrogateDataPlots(Data,models,SurrogateData,bestmodel,fitResults)

nModls = length(models);
Nsj = length(Data);

nfig=get(gcf,'Number');

mkdir figs 

%--------------------------------------------------------------------
% compare surrogate with real data 
%--------------------------------------------------------------------
nfig=nfig+1; figure(nfig);clf;

as = zeros(max([Data.Nch]),nModls,Nsj);
for sj=1:Nsj
	a(:,sj) = Data(sj).a; 

	for mdl=1:nModls;
		foo = [SurrogateData(sj).(models(mdl).name).a]; 
		%foo = reshape(foo,Data(sj).Nch,length(foo)/Data(sj).Nch);
		as(:,mdl,sj) = mean(foo==1,2);
	end

end

m = sum(a==1,2)./sum(~isnan(a),2);
s = nanstd(a==1,[],2)./sqrt(Nsj);

navat = length(unique(Data(1).avatid));
nrep = max([Data.Nch])/navat;

avatval = Data(1).avatval'*kron(eye(navat),ones(nrep,1))/nrep; 
avatval(avatval==-1)=0; avatval = avatval*2+1;
col = [1 .7 .7; .7 .7 .7; .7 1 .7];

for k=1:length(avatval);
	h=fill([.5 8.5 8.5 .5]+(k-1)*8,[0 0 1 1],col(avatval(k),:));
	hon
end
hhh = plot(sum(~isnan([Data.a])/Nsj),'k');
h(1) = plot(m,'k--','linewidth',2);
plot(find(Data(1).wordval==1),m(Data(1).wordval==1),'k.','markersize',30);
errorbar(m,s,'k.','linewidth',1);
axis tight

mas = sum(as,3)./sum(~isnan(as),3);
hh = plot(mas,'linewidth',2);


le = {models.name}; 
le = {'Data',le{:},'% responded'};
f = legend([h;hh;hhh],le,'orientation','horizontal','location','southoutside'); 
if isfield(Data(1),'word');
   set(gca,'xtick',1:length(Data(1).word),'xticklabel',Data(1).word,'xticklabelrotation',90)
end

myfig(gcf,sprintf('%s/figs/SurrogateDataPlots',fitResults));

if isfield(Data,'spin');

nfig=nfig+1; figure(nfig);clf;

	for mdl=1:nModls
		R.(models(mdl).name) = load(sprintf('%s/%s',fitResults,models(mdl).name));
		parFitEM = R.(models(mdl).name).E;
		yl = models(mdl).parnames_untr;
		[c,p] = corr(parFitEM',[Data.spin]','type','spearman'); 
		cc{mdl}=c;
		pp{mdl}=p;
		subplot(1,nModls,mdl);
			h=bar(c); 
			set(h,'facecolor',[.7 .7 .7]);
			hon
			c(p>.05)=NaN;                  h=bar(c); set(h,'facecolor',[0 1 0]);
			c(p>.05/models(mdl).npar)=NaN; h=bar(c); set(h,'facecolor',[1 0 0]);
			for k=1:models(mdl).npar
				h=text(k,min(cc{mdl}(k)+.02,.01),sprintf('p=%.3f',pp{mdl}(k)));
				set(h,'Rotation',90,'fontsize',12);
			end
			hof
			set(gca,'xtick',1:models(mdl).npar,'xticklabel',yl);
			title(models(mdl).name);
			

	end

	myfig(gcf,sprintf('%s/figs/SPINParameterCorrelations',fitResults));

end


