function labelplots(N,varargin);
%
% labelplots(N,['out'/'in'],[startlabel],[fontsize],[startsubplot])
%

if length(varargin)>=2 & ~isempty(varargin{2}); 
	k0=varargin{2};
else k0=0;
end

if length(varargin)>=3 & ~isempty(varargin{3}); 
	sz=varargin{3};
else sz=24;
end

if length(varargin)>=4 & ~isempty(varargin{4}); 
	spl0=varargin{4};
else spl0=0; 
end
%for k=(1:N)+k0
%	gtext(char(64+k+k0)),'fontsize',sz,'fontweight','bold');
%end

% get all the positions of the subplots 
t=sort(get(gcf,'children'));

% remove all but axes 
rem = [];
for k=1:length(t);
	if ~isa(t(k),'matlab.graphics.axis.Axes')
		rem = [rem k];
	end
end
t(rem) = []; 

% get positions of axes 
for k=1:N;
	pos(k,:)=get(t(k),'position');
end

% sort... 
y = unique(pos(:,2)); 
[foo,i] = sort(-pos(:,2));
pos = pos(i,:); 
t = t(i);
for k=1:length(y); 
	j = find(pos(:,2)==y(k));
	pp = pos(j,:); tt = t(j);
	[foo,i] = sort(pp(:,1)); 
	pp = pp(i,:); tt = tt(i);
	pos(j,:) = pp; 
	t(j) = tt; 
end

% for each row, find how many columns there are 
yu = unique(pos(:,2));
yn = length(yu);
for k=1:yn
	i = pos(:,2)==yu(k);
	ncol(i) = sum(i);	 % number of plots in the column 
end

if strcmp(varargin{1},'out')
	for k=spl0+1:N;
		try 
			p = t(k).OuterPosition;
			if ncol(k)==1; minfacx = .05;	
			else           minfacx = .08;
			end
			minfacy = 1.15; 
		catch
			p=pos(k,:); 
			if ncol(k)==1; minfacx = .09;	
			else           minfacx = .2;
			end
			minfacy = 1.2; 
		end
		annotation('textbox',[max(p(1)-minfacx*p(3),0) min(p(2)+minfacy*p(4),1) .0 .0],'linestyle','none','fontsize',sz,'string',char(64+k+k0));
	end 
else
	%t=sort(get(gcf,'children'));
	for k=spl0+1:N;
		%p=get(t(k),'position');
		annotation('textbox',[p(1) (p(2)+p(4)) .0 .0],'linestyle','none','fontsize',sz,'string',char(64+k+k0));
	end 

end
