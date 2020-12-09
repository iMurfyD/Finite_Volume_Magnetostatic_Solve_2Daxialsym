clear all; clc;
addpath('../solver/');

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
n = 500;
m = 500;
sdom = linspace(-8*a, 8*a, n);
zdom = linspace(-8*a, 8*a, m);
ds = sdom(2)-sdom(1);
dz = zdom(2)-zdom(1);
[XX,YY] = meshgrid(sdom,zdom);

alphaz = linspace(0.05, 0.5,30);
inf_norm = zeros(size(alphaz));

for alphi = 1:length(alphaz)

syst = struct('m',m,'n',n,'a',a,'ds',ds,'dz',dz,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,....
              'alpha', alphi);

[A,b, perm_map_debug_two_sph] = setup_system(syst);

phi_theoretical = zeros(size(XX));
u_theoretical = zeros(n*m,1);
for i = 1:length(sdom)
    for j = 1:length(zdom)
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

inf_norm(alphi) = norm(A*u_theoretical-b,Inf);

end

figure;
plot(alphaz,inf_norm, 'o-');
xlabel('\alpha');
ylabel('$||A u_t - b||_{\infty}$', 'interpreter', 'latex');
title('Error Based on Smoothed Profile Thickness');