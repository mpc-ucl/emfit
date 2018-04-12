function [D] = extractValidTrials(D)

% find out maximal depth 
dmax=0;
for sj=1:length(D);
	dmax = max(dmax,size(D(sj).choices,2));
end
Z.dmax = dmax; 

% precompute matrices for rapid inference 
% NB: assumes standard experimental setup 
Z = precomputeParams(Z);

for sj=1:length(D);

	a = D(sj).choices;
	s = D(sj).states;
	r = D(sj).rewards;

	D(sj).depths = sum(a~=0,2);
    
   d = zeros(length(r),6);
   for j=1:length(d);
       d(j,1:D(sj).depths(j)) = [D(sj).depths(j):-1:1];
   end

	% exclude no response trials 
	i = (sum(a,2)>0) & ~(r(:,1)==-200);
   D(sj).i = [1:length(a)];
   D(sj).nTotal = length(a);
   D(sj).i = D(sj).i(i);

	a = a(i,:);
	s = s(i,:);
	r = r(i,:);
	d = d(i,:);
    
	ri= r; 
	ri(ri==-70)=1; 
	ri(ri==-20)=2; 
	ri(ri== 20)=3; 
	ri(ri==140)=4; 

	an = a; an(an==0)=1;an=an-1;
	an = an*(2.^[0:dmax-1])'+1;
	sn = s(:,1);
	dn = d(:,1);

	D(sj).a  = a; % data matrix
	D(sj).s  = s; 
	D(sj).r  = r; 
	D(sj).d  = d; 
	D(sj).ri = ri; 
	D(sj).an = an; % binary code for action sequence to identify Q matrix values
	D(sj).sn = sn;
	D(sj).dn = dn;

	D(sj).Nch=length(an);
    
   % add transition matrix
   transitions = zeros(size(a));
   for ai = 1:length(a);
      transitions(ai,1:d(ai)) = s(ai,1:d(ai)) * 2 + a(ai,1:d(ai)) - 2;
   end
   D(sj).transitions = transitions;
	D(sj).Z = Z; 
end

