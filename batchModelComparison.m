function bestmodel = batchModelComparison(Data,models,resultsDir);

nModls = length(models);
Nsj = length(Data);

%--------------------------------------------------------------------
% load BICs and LLs 
%--------------------------------------------------------------------

for k=1:nModls
	try 
		loadstr = sprintf('%s/%s',resultsDir,models(k).name);
		fprintf('loading model %i %s',k,loadstr);
		R.(models(k).name) = load(loadstr);
		PL(k,:) = R.(models(k).name).stats.PL;
		LL(k,:) = R.(models(k).name).stats.LL;
		iBIC(k) = R.(models(k).name).bf.ibic;
		ilap(k) = R.(models(k).name).bf.ilap;
		fprintf('...ok\n');
	catch 
		fprintf('...ERROR, fit not loaded\n');
		iBIC(k) = NaN;
		ilap(k) = NaN;
		LL(k,:) = NaN;
	end
end

% get best model 
[foo,bestmodel] = min(iBIC);

% fitted average choice probability 
pc = exp(-LL*diag(1./[Data.Nch]));

%--------------------------------------------------------------------
% some plots 
%--------------------------------------------------------------------
fprintf('Making some plots\n');

figure(1);clf;
subplot(121)
	h=barh(mean(pc,2));set(h,'facecolor',[.7 .7 .7]);
	hon 
	plot(pc,1:nModls,'.','color',[0 0 .4]);
	plot(pc,1:nModls,'--','color',[.4 .4 1])
	hof
	set(gca,'yticklabel',{models.name},'ytick',1:nModls)
	xlabel({'Average model log posterior','probability (given prior)'})
subplot(122)
	h=barh(iBIC-min(iBIC));set(h,'facecolor',[.7 .7 .7]);
	hon
	xl = xlim; 
	plot(xl(1) + 0.8*(diff(xl)),bestmodel,'r*','markersize',20);
	for k=1:nModls;
		if k~=bestmodel
			h=text(xl(1) + .05*diff(xl),k,sprintf('\\Delta iBIC = %.1f',iBIC(k)-min(iBIC)));
			set(h,'fontsize',14)
		end
	end
	hof
		
	set(gca,'yaxislocation','right','ytick',1:nModls,'yticklabel',{models.npar});
	xlabel('\Delta iBIC score (right: #params)')
labelplots(2,'out')

mkdir([resultsDir '/figs']);
myfig(gcf,sprintf('%s/figs/ModelComparison',resultsDir));

%--------------------------------------------------------------------
% output as latex table 
%--------------------------------------------------------------------
fprintf('Outputting model comparison data as latex table\n');
fid = fopen(sprintf('%s/modelComparison.tex',resultsDir),'w');
fprintf(fid,'model name & \\# params & iBIC & Description');
for k=1:nModls
	try 
		fprintf(fid,'\n\\\\%s & %i & %g & %s',models(k).name,models(k).npar,R.(models(k).name).bf.ibic,models(k).descr);
	catch 
		fprintf(fid,'\n\\\\%s & %i &    & %s',models(k).name,models(k).npar,models(k).descr);
	end
end
fclose(fid);

