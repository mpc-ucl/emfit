function p = wfpt_dz_high_boundary(t,v,a,z,err)
% This function calculates the probability for Wiener first passage time
% according to Navarro et al 2009.
% t = time [s]
% v = drift
% a = boundary
% z = starting point
% err = 10^-29 (default)
if t<0
    p=(10^-20)^(-t);
    return      
end

tt = t/(a^2); % use normalized time
w = z/a; % convert to relative start point --> since sp, i.e. z = spf*a,
% can ignore the boundary here

% calculate number of terms needed for large t
if pi*tt*err < 1 % if error threshold is set low enough
    kl = sqrt(-2*log(pi*tt*err)./(pi^2*tt)); % bound
    kl = max(kl,1/(pi*sqrt(tt))); % ensure boundary conditions met
else % if error threshold set too high
    kl = 1/(pi*sqrt(tt)); % set to boundary condition
end

% calculate number of terms needed for small t
if 2*sqrt(2*pi*tt)*err<1 % if error threshold is set low enough
    ks = 2+sqrt(-2*tt.*log(2*sqrt(2*pi*tt)*err)); % bound
    ks = max(ks,sqrt(tt)+1); % ensure boundary conditions are met
else % if error threshold was set too high
    ks=2; % minimal kappa for that case
end

% compute f(tt|0,1,w)
p2 = 0; %initialize density
p2_dz = 0; 
if ks < kl % if small t is better...
    K = ceil(ks); % round to smallest integer meeting error
    for k = -floor((K-1)/2):ceil((K-1)/2) % loop over k
        p2 = p2 + (w + 2*k)*exp(-((w + 2*k)^2)/2/tt); % increment sum
        p2_dz = p2_dz + w*exp(-((w + 2*k)^2)/2/tt) + (w + 2*k)*(-2/2/tt)*(w + 2*k)*w*exp(-((w + 2*k)^2)/2/tt);
    end
    p2 = p2/sqrt(2*pi*tt^3); % add constant term
else % if large t is better...
    K=ceil(kl); % round to smallest integer meeting error
    for k = 1:K
        p2 = p2+k*exp(-(k^2)*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
        p2_dz = p2_dz+k*exp(-(k^2)*(pi^2)*tt/2)*(k*pi*w)*cos(k*pi*w);
    end
    p2 = p2*pi; % add constant term
    p2_dz = p2_dz*2;
end

% convert to f(t|v,a,w)
p1 = exp(-v*a*w -(v^2)*t/2)/(a^2);
p1_dz = (-v*a*w)*exp(-v*a*w -(v^2)*t/2)/(a^2);

% make complete derivative
p = p1*p2_dz+p1_dz*p2; 

if p < 10^-20
    p = 10^-20*tt;
end

end
