function [name, file]=myfig(gcf,name,mfile,varargin);
% [name, figurefile]=myfig(gcf,name,mfile,filelist);
% take input gcf and set paperpositionmode to auto and all tickmodes to manual
% name should include directory 
% print an eps file and a fig file
% additionally copy the m.file used to make the figure there
% and finally also copy all files in the filelist cell to a .tar.gz in the same
% place 

in=input(['do you really want to save the figure as ',name,'.eps? [y for yes, else no]'],'s');

if length(varargin)==1
	furtherfiles = varargin{1};
end

if strcmp(in,'y')

	set(gcf,'paperpositionmode','auto')
	t=get(gcf,'children');
	for k=1:length(t);
		try 
			set(t(k),'xtickmode','manual','ytickmode','manual','ztickmode','manual');
		end
	end

	print(gcf,'-depsc2',strcat(name,'.eps'))
	saveas(gcf,strcat(name,'.fig'));

	if exist('mfile'); eval(['!cp ',mfile,' ',name,'.m']);end
	if exist('furtherfiles'); 
		filelist=[name '.eps ' name '.fig ' file ' '];
		for k=1:length(furtherfiles)
			filelist=[filelist ' ' furtherfiles{k}];
		end
		eval(['! tar -zcvf ' name '.tar.gz ' filelist]);
	end

end
