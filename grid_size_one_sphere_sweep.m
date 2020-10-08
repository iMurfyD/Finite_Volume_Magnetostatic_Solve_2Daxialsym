clear all; clc;

%% Set up parameters
a = 1.4e-6; % Radius of sphere
sep = 2.2*a;
r1 = [0 sep/2]';
r2 = [0 sep/2]'; 
perm_free_space = 1;%4*pi*1.00000000082e-7; % H*m^-1 
                                         % permiability of free space
                                         
susc = 0.96  ; % Suscepibility of material, arbitrary

perm = perm_free_space*(1+susc); % Linear media, eqn 6.30 Griffiths
H0 = [0 477]'; % A/m
%% Set up the grid (n x m 2D grid)
gsizes = [200 300 400 500 600 700 800 900 1000 1100 1200 1300 1500 ...
          1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 ...
          2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 ...
          4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5000]';


alpha = 0.27;
inf_norm = zeros(size(gsizes));

for gridi = 1:length(gsizes)
n = gsizes(gridi);
m = gsizes(gridi);
xdom = linspace(-16*a, 16*a, n);
ydom = linspace(-16*a, 16*a, m);

dx = xdom(2)-xdom(1);
dy = ydom(2)-ydom(1);
[XX,YY] = meshgrid(xdom,ydom);

syst = struct('n',n,'m',m,'a',a,'dx',dx,'dy',dy,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,...
              'alpha', alpha);

[A,b, perm_map_debug_two_sph] = setup_system_sparse(syst);

phi_theoretical = zeros(size(XX));
u_theoretical = zeros(n*m,1);
for i = 1:length(xdom)
    for j = 1:length(ydom)
        idx = n*(i-1) + j; % Ordering of nodes in 1D
        rvec = [XX(i,j)-r1(1) YY(i,j)-r1(2)]';
        r = norm(rvec,2);
        thta = atan2(rvec(1),rvec(2));
        if(r >= a) % Outside
        phi_theoretical(i,j) = -H0(2)*r*cos(thta) + ...
                           H0(2)*power(a,3)*((perm-perm_free_space)/...
                                          (perm+2*perm_free_space))*...
                                        cos(thta)/power(r,2);
        else % Inside
            phi_theoretical(i,j) = -H0(2)*((3*perm_free_space)/...
                                           (perm+2*perm_free_space))*...
                                     r*cos(thta);
    
        end
        u_theoretical(idx) = phi_theoretical(i,j);
    end
end

inf_norm(gridi) = norm(A*u_theoretical-b,Inf);

end

figure;
plot(gsizes,inf_norm, 'o-');
xlabel('n x n Grid Size');
ylabel('$||A u_t - b||_{\infty}$', 'interpreter', 'latex');
title('Error Based on Smoothed Profile Thickness');