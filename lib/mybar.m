function h=mybar(varargin)
%
% function h=mybar(x,m,c)
%
% make bar, and change face color
% either two or three arguments: 
% if two arguments: m and c
% if three arguments: additionally x
% 
%  c can have the following dimensions: 
%
%    1x1: all bars same gray c
% 	  1x3: all bars same colour [c(1) c(2) c(3)]
% 	  nx1: each bar different gray c(k)
% 	  nx3: each bar different colour [c(k,1) c(k,2) c(k,3)]
%
%	Quentin Huys 2010
%  qhuys@cantab.net

	if     nargin==2;                  m = varargin{1}; c = varargin{2};
	elseif nargin==3; x = varargin{1}; m = varargin{2}; c = varargin{3};
	else   error('not the right  number of inputs')
	end


if min(size(m))==1 & max(size(m))==size(c,1)	 % each bar with different colour 
	
	nb = max(size(m));
	if nargin == 2; x = 1:nb; end
	if min(size(c))==1; c = c(:)*ones(1,3);end

	if ishold; holdmem=1; else; holdmem=0; cla; end

	hold on	

	for k=1:nb
		h(k) = bar(x(k),m(k));
		set(h(k),'facecolor',c(k,:));
	end

	if holdmem == 0; hold off;end

else	% all bars of a group with same colour 

	if     nargin==2; h=bar(m);
	elseif nargin==3; h=bar(x,m);
	end

	if     (size(c,1)==1 & size(c,2)==1); set(h,'facecolor',c*ones(1,3));
	elseif (size(c,1)==1 & size(c,2)==3); set(h,'facecolor',c);
	elseif (size(c,1)>1  & size(c,2)==3); for k=1:size(c,1); set(h(k),'facecolor',c(k,:));         end
	elseif (size(c,1)>1  & size(c,2)==1); for k=1:size(c,1); set(h(k),'facecolor',c(k)*ones(1,3)); end
	else	 error('color matrix size is wrong')
	end

end
