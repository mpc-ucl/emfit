% plot split-half reliability for emfit models

% Note: only does this for the first model at the moment

i = 1;
try
   file1ToLoad = sprintf('fitResults/%s_firstHalf.mat', models(i).name); 
   load(file1ToLoad)
    
   E_fh = E;
   
   file2ToLoad = sprintf('fitResults/%s_secondHalf.mat', models(i).name); 
   load(file2ToLoad)
   
   E_sh = E;

   Np = models(i).npar;
   figure;
   x1 = 1;
   pi=1:size(E_fh,2);
   for pr = 1:Np
        subplot(x1,Np,pr); 
		plot(E_fh(pr,:), E_sh(pr,:),'k.','markersize',20);
		hon
		plot(E_fh(pr,~pi),E_fh(pr,~pi),'r.','markersize',20);
		yl=ylim; xl=xlim; plot(xl,[0 0],'k'); plot([0 0],yl,'k');
		plot([-10 10],[-10 10],'k');
		ylim(yl);xlim(xl);
		hof
		[c,p] = corr(E_fh(pr,pi)',E_sh(pr,pi)','type','spearman');
		title({sprintf('Model: %s', models(i).name), sprintf('Split-half correlation parameter %i',pr),sprintf('amongst ok fit: p=%.1g r=%.2g',p,c)});
		xlabel('1st half');
		ylabel('2nd half');
    end
  
   
catch
    fprintf('not all data saved to compute split-half reliability')    
end



try
   file1ToLoad = sprintf('fitResults/%s_odd.mat', models(i).name); 
   load(file1ToLoad)
    
   E_fh = E;
   
   file2ToLoad = sprintf('fitResults/%s_even.mat', models(i).name); 
   load(file2ToLoad)
   
   E_sh = E;

   Np = models(i).npar;
   figure;
   x1 = 1;
   pi=1:size(E_fh,2);
   for pr = 1:Np
        subplot(x1,Np,pr); 
		plot(E_fh(pr,:), E_sh(pr,:),'k.','markersize',20);
		hon
		plot(E_fh(pr,~pi),E_fh(pr,~pi),'r.','markersize',20);
		yl=ylim; xl=xlim; plot(xl,[0 0],'k'); plot([0 0],yl,'k');
		plot([-10 10],[-10 10],'k');
		ylim(yl);xlim(xl);
		hof
		[c,p] = corr(E_fh(pr,pi)',E_sh(pr,pi)','type','spearman');
		title({sprintf('Model: %s', models(i).name),sprintf('Split-half correlation parameter %i',pr),sprintf('amongst ok fit: p=%.1g r=%.2g',p,c)});
		xlabel('odd');
		ylabel('even');
    end
  
   
catch
    fprintf('not all data saved to compute split-half reliability\n')    
end