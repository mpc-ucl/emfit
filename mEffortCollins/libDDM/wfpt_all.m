function [p, dv, da, dz, dt, dv_ps, dz_ps, da_ps, dt_ps]  = wfpt_all(t,v,a,z,err)
% This function calculates the probability for Wiener first passage time
% according to Navarro et al 2009.
% t = time 
% v = drift
% a = boundary
% z = starting point
% err = 10^-29 (default)
% if t<0
%     p=(10^-20)^(-t);
%     dz=0; 
%     da=0; 
%     dv=0; 
%     dt=log(10^20); 
%     dv_ps = 0; 
%     dz_ps = 0; 
%     da_ps = 0; 
%     dt_ps = 4^(10*t+1)*5^(20*t+1)*log(10); 
%     return      
% end

tt = t/(a^2); % use normalized time
w = z/a; % convert to relative start point

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
p = 0; %initialize density
dz = 0; 
da = 0; 
dt = 0; 
if ks < kl % if small t is better...
    K = ceil(ks); % round to smallest integer meeting error
    for k = -floor((K-1)/2):ceil((K-1)/2) % loop over k
%      if ((z/a) + 2*k)^2/(2*(t/(a^2)))<0.00135
%          p=0;
%          dv=-z-v*t;
%          da=0;
%          dz=-1e4;
%          dt=1e4;
%          return      
%      end

        p = p + ((z/a) + 2*k)*exp(-(((z/a) + 2*k)^2)/(2*(t/(a^2)))); % increment sum 
        dz = dz+((1/a)*exp(-(((z/a) + 2*k)^2)/(2*(t/(a^2))))+((z/a) + 2*k)*exp(-(((z/a) + 2*k)^2)/(2*(t/(a^2))))*(-1/(t/a^2))*((z/a) + 2*k)*(1/a));   
        da = da+(-1)*exp(-(2*k*a+z)^2/(2*t))*(8*k^3*a^3+8*k^2*a^2*z+2*k*a*z^2+t*z)/(t*a^2); 
        dt = dt+a^2*((z/a) + 2*k)^3*exp(-(((z/a) + 2*k)^2)/(2*(t/(a^2))))/(2*t^2);
    end
    dz_part1 = dz/sqrt(2*pi*(t/(a^2))^3);
    da_part11 = da;
    dt_part11 = dt; 
    dz = dz * (1/p);
    da = da*(1/p); 
    dt = dt*(1/p); 
    p_part11 = p;
    p = p/sqrt(2*pi*(t/(a^2))^3); % add constant term 
    p_part12 = 1/sqrt(2*pi*(t/(a^2))^3);
    p_part1 = p;
    dv_part1 = p;
    da = da+(3/a); 
    dt = dt-(3/(2*t)); 
    da_part12 = 3*t^3/(sqrt(2*pi)*a^7*(t^3/a^6)^(3/2));
    dt_part12 = -3*t^2/(2*sqrt(2*pi)*a^6*(t^3/a^6)^(3/2));
    da_part1 = da_part11*p_part12+da_part12*p_part11;
    dt_part1 = dt_part11*p_part12+dt_part12*p_part11; 
    
else % if large t is better...
    K=ceil(kl); % round to smallest integer meeting error
    for k = 1:K
        p = p+k*exp(-(k^2)*(pi^2)*(t/(a^2))/2)*sin(k*pi*(z/a)); % increment sum
        dz = dz + k*exp(-(k^2)*(pi^2)*(t/(a^2))/2)*cos(k*pi*(z/a))*(k*pi/a); 
        da = da+k*(exp(-(k^2)*(pi^2)*(t/(a^2))/2)*((k^2)*(pi^2)*(t/(a^3)))*sin(k*pi*(z/a))+exp(-(k^2)*(pi^2)*(t/(a^2))/2)*cos(k*pi*(z/a))*(-k*pi*z/a^2));
        dt = dt+k*exp(-(k^2)*(pi^2)*(t/(a^2))/2)*sin(k*pi*(z/a))*(-(k^2)*(pi^2)*(1/(2*a^2)));
    end
    dz_part1 = dz*pi; 
    da_part1 = da; 
    dt_part1 = dt; 
    dz = dz*(1/p); 
    da = da*(1/p); 
    dt = dt*(1/p); 
    p = p*pi; % add constant term
    p_part1 = p; 
    dv_part1 = p; 
    da_part1 = da_part1*pi; 
    dt_part1 = dt_part1*pi; 
end

% convert to f(t|v,a,w)
p = p*exp(-v*z-(v^2)*t/2)/(a^2);
dv = -z-v*t;
dv_ps = (-z-v*t)*dv_part1*exp(-v*z-(v^2)*t/2)/(a^2);
da = (-2/a)+da;
da_ps = exp(-v*z-(v^2)*t/2)*(-2)*a^(-3)*p_part1+da_part1*exp(-v*z-(v^2)*t/2)/(a^2);
dz = -v+dz; 
dz_ps = (exp(-v*z-(v^2)*t/2)/(a^2))*dz_part1+p_part1*(exp(-v*z-(v^2)*t/2)/(a^2))*(-v); 
dt = (-v^2/2)+dt; 
dt_ps = (exp(-v*z-(v^2)*t/2)/(a^2))*(-v^2/2)*p_part1+dt_part1*(exp(-v*z-(v^2)*t/2)/(a^2));


% if p < 10^-20
%     p = 10^-20*tt;
%     dv = 0; 
%     dz = 0; 
%     da = 10^-20*t*(-2)/a^3;
%     dt = 0; 
% end

end
