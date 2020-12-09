%% Exact sphere method, finite volume
% Ian DesJardin
% August 2020

clear all; close all; clc;
addpath('../solver/');

%% Set up parameters
a = 1.4e-6; % Radius of sphere
sep = 2.2*a;
r1 = [0 0]';
r2 = [0 0]'; 
perm_free_space = 4*pi*1.00000000082e-7; % H*m^-1 
                                         % permiability of free space
                                         
susc = 0.96; % Suscepibility of material, arbitrary

perm = perm_free_space*(1+susc); % Linear media, eqn 6.30 Griffiths
H0 = [0.0 477]'; % A/m

%% Define what stuff to plot
four_plots = 0;
interp_match_du_paper = 1   ;
perm_map_debug=1;
spy_mat = 0;
hy_trifold = 1;
hx_trifold = 0;
hx_trifold_no_normalization =0;
hy_trifold_no_normalization=0;
phi_trifold=0;
hmag=1;
debug_force_calc = 1;  
one_sph_debug = 1;

%% Set up the grid (n x m 2D grid)
n = 1200;
m = 1200;
sdom = linspace(-8*a, 8*a, n);
zdom = linspace(-8*a, 8*a, m);
ds = sdom(2)-sdom(1);
dz = zdom(2)-zdom(1);
[XX,YY] = meshgrid(sdom,zdom);

syst = struct('m',m,'n',n,'a',a,'ds',ds,'dz',dz,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,....
              'alpha', 0.2517);

%% Form FV Matrix
fprintf('Setting up Linear System\n');
[A,b, perm_map_debug_two_sph] = setup_system(syst);
fprintf('Solving Linear System with \\ Operator\n');
u = A\(b);
fprintf('Done Solving Linear System with \\ Operator\n');

phi = spread_1D_into_2D(u,syst);

if(spy_mat)
figure;
spy([A b]);
title('Two Sphere [A b] Matrix');

figure;
pcolor(A);
title('pcolor of A');
colorbar;
set(gca, 'YDir','reverse')
end

if(one_sph_debug)
% Find theoretical 1 sphere phi
phi_theoretical = zeros(size(XX));
u_theoretical = zeros(n*m,1);
for i = 1:length(sdom)
    for j = 1:length(zdom) % Because Meshgrid index y, x, z 
        idx = indexf(i,j,syst);
        rvec = [XX(j,i)-r1(1) YY(j,i)-r1(2)]';
        r = norm(rvec,2);
%         thta = acos(rvec(3)/r);
        cthta = rvec(2)/r;
        if(r >= a) % Outside
        pt = (perm-perm_free_space)/(perm+2*perm_free_space);
        phi_theoretical(i,j) = -H0(2)*r*cthta + ...
                           H0(2)*power(a,3)*pt*cthta/power(r,2);
        else % Inside
            phi_theoretical(i,j) = -H0(2)*((3*perm_free_space)/...
                                           (perm+2*perm_free_space))*...
                                     r*cthta;
    
        end
        u_theoretical(idx) = phi_theoretical(i,j);
    end
end

[HXt, HYt] = gradient(phi_theoretical,sdom,zdom);
tmp = HXt;
HXt=-HYt;
HYt=-tmp';

[HXn,HYn] = gradient(phi, sdom,zdom);    
HXn=-HXn;
HYn=-HYn;

%% Plotting
figure; subplot(1,3,1);
absHn = sqrt(HXn.^2+HYn.^2);
absHt = sqrt(HXt.^2+HYt.^2);
pc = pcolor(XX./a,YY./a,log10(abs(absHn-absHt))); set(pc, 'EdgeColor', 'none');
colorbar; title('log_{10}(|\Delta H|)');
axis equal;
xlim([-2 2]); ylim([-2 2]);

subplot(1,3,2);
pc = pcolor(XX./a,YY./a,log10(abs(HXn-HXt))); set(pc, 'EdgeColor', 'none');
colorbar; title('log_{10}(|\Delta H_{x}|)');
axis equal;
xlim([-2 2]); ylim([-2 2]);

subplot(1,3,3);
pc = pcolor(XX./a,YY./a, log10(abs(HYn-HYt))); 
set(pc, 'EdgeColor', 'none');
colorbar; title('log_{10}(|\Delta H_{y}|)');
axis equal;
xlim([-2 2]); ylim([-2 2]);

figure;subplot(1,3,1);
pc = pcolor(XX./a,YY./a,absHn); set(pc, 'EdgeColor', 'none');
colorbar; title('|H| Calculated');
axis equal;
xlim([-2 2]); ylim([-2 2]);
subplot(1,3,2);
pc = pcolor(XX./a,YY./a,HXn); set(pc, 'EdgeColor', 'none');
colorbar; title('H_{x} Calculated');
axis equal;
xlim([-2 2]); ylim([-2 2]);
subplot(1,3,3);
pc = pcolor(XX./a,YY./a,HYn); set(pc, 'EdgeColor', 'none');
colorbar; title('H_{y} Calculated');
axis equal;
xlim([-2 2]); ylim([-2 2]);

figure; subplot(1,3,1);
absHn = sqrt(HXn.^2+HYn.^2);
absHt = sqrt(HXt.^2+HYt.^2);
pc = pcolor(XX./a,YY./a,absHn-absHt); set(pc, 'EdgeColor', 'none');
colorbar; title('\Delta H');
axis equal;
xlim([-2 2]); ylim([-2 2]);

subplot(1,3,2);
pc = pcolor(XX./a,YY./a,HXn-HXt); set(pc, 'EdgeColor', 'none');
colorbar; title('\Delta H_{x}');
axis equal;
xlim([-2 2]); ylim([-2 2]);

subplot(1,3,3);
pc = pcolor(XX./a,YY./a, HYn-HYt); 
set(pc, 'EdgeColor', 'none');
colorbar; title('\Delta H_{y}');
axis equal;
xlim([-2 2]); ylim([-2 2]);

end    
