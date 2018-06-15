function h=mydeb(xx,y,v,varargin);
%
% MYDEB  Plot double errorbars, with 1 st err in red and 95% conf interval in green. 
%
% H=MYDEB(X,Y,V,...);
% 
% X is X to plot against. If X=0, then generates appropriate X. 
% Y is individual means
% V is individual standard deviations or standard errors
%
% 'funct' is optional. If given, then the data in Y and V are first passed
% through the function, e.g. a sigmoid. 
%
% Quentin Huys 2010. 
% qhuys@cantab.net

DX=[.22 .15 .22 .275 .31 linspace(.332,.367,7)];

if xx==0
	[nx,ny] = size(y);
	if ~any([nx ny]==1)
		dx = DX(ny);
		xx = (1:nx)'*ones(1,ny) + ones(nx,1)*linspace(-dx,dx,ny);
	else
		xx = 1:max(nx,ny);
	end
end

if ishold; holdmem=1; else; holdmem=0;end

hold on
ci=1; se=1;
nv = length(varargin);
if nv>0
	for k=1:2:nv
		switch varargin{k}
			case 'funct'
				fun=varargin{k+1};
			case 'cionly';
				se=0;
				fun=inline('x');
			case 'seonly';
				ci=0;
				fun=inline('x');
		end
	end

	m = fun(y);
	c1 = fun(y-1.96*v)-m;
	c2 = fun(y+1.96*v)-m;
	s1 = fun(y-v)-m;
	s2 = fun(y+v)-m;

	if ci; errorbar(xx,m,c1,c2,'g','linewidth',2,'linestyle','none');end
	if se; errorbar(xx,m,s1,s2,'r','linewidth',2,'linestyle','none');end

else
	if ci; errorbar(xx,y,1.96*v,'g','linewidth',2,'linestyle','none');end
	if se; errorbar(xx,y,     v,'r','linewidth',2,'linestyle','none');end
end

if holdmem == 0; hold off;end
