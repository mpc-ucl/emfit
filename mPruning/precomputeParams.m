function [Z] = precomputeParams(Z);

fprintf('standard reward/transition matrix setup assumed!\n');

R_matrix = [140 20; -20 -70; -70 -20; 20 -20; -70 -20; 20 -20];
T_matrix = [2  4;   3   5;   6   4;  2   5;   1   6;  3   1]; 
trans_matrix = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12]; 

Z.Rm = R_matrix;
Z.Tm = T_matrix;
Z.r3=min(Z.Rm(:));

tmp = unique(Z.Rm(:));
Z.Rmi = Z.Rm;
for k=1:4; 
	Z.Rmi(Z.Rm==tmp(k)) = k; 
	tmp2 = 1:4; tmp2(k)=[];
	Z.Rmo(k,:) = tmp2;
end
clear tmp tmp2

% in the following matrices max dep 5 was hardcoded
dmax = Z.dmax;
%Z.Hd = Inf*ones(2^dmax,2^dmax,dmax-2); 
Z.Ps = zeros(dmax,2^dmax,dmax-2,6);
Z.Pg = zeros(dmax,2^dmax,dmax-2,6);
Z.Psgd = zeros(dmax,2^dmax,dmax-2,6);
Z.Pggd = zeros(dmax,2^dmax,dmax-2,6);
for d=3:dmax % max dep 5 was hardcoded
	% get all possible action sequences 
	clear tmp 
	tmp = dec2bin(0:2^d-1); 
	for k=1:d; Z.A(d-k+1,1:2^d,d-2)=str2num(tmp(:,k))+1;end

	% for each state, apply all action sequences and get state and reward
	% sequences 
	for s=1:6
		for k=1:2^d
			aa=Z.A(1:d,k,d-2);
			ss=s;
			Z.S(1,k,d-2,s)=s;
            aidoffset = 2;
			for dd=1:d
				r  = Z.Rm (ss,aa(dd));
				ri = Z.Rmi(ss,aa(dd));
				ss = Z.Tm (ss,aa(dd));
				Z.R (dd  ,k,d-2,s) = r; % reward matrix for all state-action-seqs
				Z.Ri(dd  ,k,d-2,s) = ri; % reward identity matrix
                % transition id:
				Z.trans (dd  ,k,d-2,s) = trans_matrix (ss,aa(dd));
                
                aidx = aa(dd) - 1 + aidoffset;
                Z.actionIdx (dd,k,d-2,s) = aidx;
                aidoffset = aidx*2;
                
                
				for kk=1:4; Z.Rikk(dd,k,kk,d-2,s) = ri==kk; end
				Z.S (dd+1,k,d-2,s) = ss; % states
				Z.D (dd  ,k,d-2,s) = dd-1; 
				if r==Z.r3; Z.Ps(dd+1:d,k,d-2,s) = Z.Ps(dd+1:d,k,d-2,s)+1; % exponents of specific gamma
				else      ; Z.Pg(dd+1:d,k,d-2,s) = Z.Pg(dd+1:d,k,d-2,s)+1; % exponents of general gamma
				end
				if r==Z.r3
					if dd==1 
						Z.Psgd(dd+1:d,k,d-2,s) = Z.Psgd(dd+1:d,k,d-2,s)+1;
					elseif dd>1 & Z.R(dd-1,k,d-2,s)~=Z.r3
						Z.Psgd(dd+1:d,k,d-2,s) = Z.Psgd(dd+1:d,k,d-2,s)+1;
					else 
						Z.Pggd(dd+1:d,k,d-2,s) = Z.Pggd(dd+1:d,k,d-2,s)+1;
					end
				else
					Z.Pggd(dd+1:d,k,d-2,s) = Z.Pggd(dd+1:d,k,d-2,s)+1;
                end
			end
		end
	end

%	% get Hamming distances between all action pairs
%	a = Z.A(1:d,1:2^d,d-2)==1;
%	for k=1:2^d
%        for j=1:2^d
%            Z.Hd(k,j,d-2) = sum(xor(a(:,k),a(:,j)));
%        end
%    end
    
end


