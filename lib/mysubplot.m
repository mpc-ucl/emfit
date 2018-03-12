function ax=mysubplot(J,K,varargin)

% mysubplot(y,x,x0,y0,xl,yl);
if exist('varargin'); v=varargin;end

if 	 length(varargin)==0; yl=.85; xl=.85; y0=.1;  x0=.1;  
elseif length(varargin)==1; yl=.85; xl=.85; y0=.1;  x0=v{1};
elseif length(varargin)==2; yl=.85; xl=.85; y0=v{2};x0=v{1};
elseif length(varargin)==3; yl=.85; xl=v{3};y0=v{2};x0=v{1};
elseif length(varargin)>=4; yl=v{4};xl=v{3};y0=v{2};x0=v{1};
end
if length(varargin)>=5; xspacing=v{5};else; xspacing=0.1; end
if length(varargin)>=6; yspacing=v{6};else; yspacing=0.1; end

%clf;
l=0;
for j=1:J;
	for k=1:K;l=l+1;
		ax(l)=axes('position',[x0+(k-1)*xl/K y0+(J-j)*yl/J xl/(K+xspacing) yl/(J+yspacing)]);
		if k>1; set(gca,'yticklabel',[]);end
		if j~=J; set(gca,'xticklabel',[]);end
	end
end
