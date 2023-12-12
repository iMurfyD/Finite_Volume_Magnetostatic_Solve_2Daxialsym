%% Exact sphere method, finite volume
% Ian DesJardin
% August 2020

clear all; close all; clc;
addpath('../solver/');
 
%% Set up parameters
a = 1.0e-2; % Radius of sphere, meters
sep = 2.005*a;
r1 = [0 -sep/2]';
r2 = [0 sep/2]'; 
perm_free_space = 4*pi*1.00000000082e-7; % H*m^-1 
                                         % permiability of free space
                                         
susc = 5000; % Suscepibility of material, arbitrary

perm = perm_free_space*(1+susc); % Linear media, eqn 6.30 Griffiths
H0 = [0.0 2000e-9/perm_free_space]'; % A/m

%% Define what stuff to plot
four_plots = 0;
interp_match_du_paper = 0  ;
perm_map_debug=1;
spy_mat = 0;
hy_trifold = 0;
hx_trifold = 0;
hx_trifold_no_normalization =0;
hy_trifold_no_normalization=0;
phi_trifold=0;
hmag=0;
debug_force_calc = 0;  
one_sph_debug = 0;
for_paper = 1;

%% Set up the grid (n (x,s) x m (y,z) 2D axial grid)
n = 6000;
m = 6000;
sdom = linspace(-8*a, 8*a, n);
zdom = linspace(-8*a, 8*a, m);
ds = sdom(2)-sdom(1);
dz = zdom(2)-zdom(1);
[XX,YY] = meshgrid(sdom,zdom);

syst = struct('n',n,'m',m,'a',a,'ds',ds,'dz',dz,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,...
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

[HXn,HYn] = gradient(phi, sdom, zdom);    
HXn=-HXn;
HYn=-HYn;

%% Plotting

figure; subplot(1,3,1);
absHn = sqrt(HXn.^2+HYn.^2);
pc = pcolor(XX./a,YY./a,absHn); set(pc, 'EdgeColor', 'none');
colorbar; title('|H|');colormap hot;
axis equal;
xlim([-3 3]); ylim([-3 3]);

subplot(1,3,2);
pc = pcolor(XX./a,YY./a,HXn);
set(pc, 'EdgeColor', 'none');
colorbar; title('H_{x}');colormap hot;
axis equal;
xlim([-3 3]); ylim([-3 3]);

subplot(1,3,3);
pc = pcolor(XX./a,YY./a,HYn); set(pc, 'EdgeColor', 'none');
colorbar; title('H_{y}'); colormap hot;
axis equal;
xlim([-3 3]); ylim([-3 3]);


%% Form FV Matrix and solve for just left sphere
r2 = r1;

syst = struct('n',n,'m',m,'a',a,'ds',ds,'dz',dz,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,...
              'alpha', 0.2517);
          
[A,b, perm_map_debug_left_sph] = setup_system(syst);
if(spy_mat)
figure;
spy([A b]);
title('Left Sphere A Matrix');
end
u = A\b;
phi_left = spread_1D_into_2D(u,syst);
perm_map_debug_left_sph = spread_1D_into_2D(perm_map_debug_left_sph, syst);
[HX_left, HY_left] = gradient(phi_left,sdom,zdom);
HX_left = -HX_left;
HY_left = -HY_left;

%% Form FV Matrix and solve for just right sphere
r2 = [0 sep/2]';
r1 = r2;

syst = struct('n',n,'m',m,'a',a,'ds',ds,'dz',dz,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,...
              'alpha', 0.2517);
          
[A,b, perm_map_debug_right_sph] = setup_system(syst);
if(spy_mat)
figure;
spy([A b]);
title('Right Sphere A Matrix');
end
u = A\b;
phi_right = spread_1D_into_2D(u,syst);
perm_map_debug_right_sph = spread_1D_into_2D(perm_map_debug_right_sph, syst);
[HX_right, HY_right] = gradient(phi_right,sdom,zdom);
HX_right=-HX_right;
HY_right=-HY_right;


%% Formulate the maxwell stress tensor and integrate force around sphere
%  Not the recommended way to do this numerically
f1 = [0 0]';
f2 = [0 0]';
perm_map_debug_two_sph = spread_1D_into_2D(perm_map_debug_two_sph,syst);
[FMX, FMY] = gradient(perm_map_debug_two_sph,sdom,zdom);
FMX = -0.5*(HXn.^2+HYn.^2).*FMX;
FMY = -0.5*(HXn.^2+HYn.^2).*FMY;

if(debug_force_calc)    
figure; 
pc = pcolor(XX,YY,FMX); set(pc, 'EdgeColor', 'none');
xlim([-3e-6 3e-6]); ylim([-3e-6 3e-6]);
title('FMX');
colorbar; 

figure; 
pc = pcolor(XX,YY,FMY); set(pc, 'EdgeColor', 'none');
hold on; % Don't let plot blow away the image.
theta = 0 : 0.01 : 2*pi;
radius = a;
xx = radius * cos(theta) + r1(1);
yy = radius * sin(theta) + r1(2);
plot(xx, yy, 'r-', 'LineWidth', 1);
theta = 0 : 0.01 : 2*pi;
radius = a;
xx = radius * cos(theta) + r2(1);
yy = radius * sin(theta) + r2(2);
plot(xx, yy, 'r-', 'LineWidth', 1);
xlim([-3e-6 3e-6]); ylim([-3e-6 3e-6]);
title('FMY');
colorbar;
end

c_map = zeros(size(XX));
fx_per_vol_dV = zeros(size(XX));
fy_per_vol_dV = zeros(size(XX));
test_fV = 1/((4/3)*pi*a^3);
for ii = 1:m
    for jj = 1:n
        if(YY(ii,jj) <= 0) % r1
            if(abs(XX(ii,jj)) <= mean([ds dz]))
                circum = 2*pi*mean([ds dz])/2;
            else
                circum = 2*pi*abs(XX(ii,jj))/2;
            end
            c_map(ii,jj) = circum;
            f_per_vol_dV = [FMX(ii,jj) FMY(ii,jj)]'*ds*dz*abs(circum);
            fx_per_vol_dV(ii,jj) = f_per_vol_dV(1);
            fy_per_vol_dV(ii,jj) = f_per_vol_dV(2);  
            f1 = f1 + f_per_vol_dV;
%             end
        else % r2
            if(abs(XX(ii,jj)) <= mean([ds dz]))
                circum = mean([ds dz]);
            else
                circum = 2*pi*abs(XX(ii,jj))/2;
            end
            c_map(ii,jj) = circum;
            f_per_vol_dV = [FMX(ii,jj) FMY(ii,jj)]'*ds*dz*abs(circum);
            fx_per_vol_dV(ii,jj) = f_per_vol_dV(1);
            fy_per_vol_dV(ii,jj) = f_per_vol_dV(2);  
            f2 = f2 + f_per_vol_dV;
%             end
        end
       
    end
end

if(debug_force_calc)
    
fprintf('mean c = %g\n', mean(circum,'all'));
fprintf('mean c/(2*pi) = %g\n', mean(circum,'all')/(2*pi));


figure; 
pc = pcolor(XX,YY,c_map); set(pc, 'EdgeColor', 'none');
hold on; % Don't let plot blow away the image.
theta = 0 : 0.01 : 2*pi;
radius = a;
xx = radius * cos(theta) + r1(1);
yy = radius * sin(theta) + r1(2);
plot(xx, yy, 'r-', 'LineWidth', 1);
xlim([-3e-6 3e-6]); ylim([-3e-6 3e-6]);
theta = 0 : 0.01 : 2*pi;
radius = a;
xx = radius * cos(theta) + r2(1);
yy = radius * sin(theta) + r2(2);
plot(xx, yy, 'r-', 'LineWidth', 1);
xlim([-3e-6 3e-6]); ylim([-3e-6 3e-6]);
title('depth');
colorbar; 
caxis([0 pi*a]);

figure; 
pc = pcolor(XX,YY,fx_per_vol_dV); set(pc, 'EdgeColor', 'none');
xlim([-3e-6 3e-6]); ylim([-3e-6 3e-6]);
title('fx_per_vol');
colorbar; 

figure; 
pc = pcolor(XX,YY,fy_per_vol_dV); set(pc, 'EdgeColor', 'none');
hold on; % Don't let plot blow away the image.
theta = 0 : 0.01 : 2*pi;
radius = a;
xx = radius * cos(theta) + r1(1);
yy = radius * sin(theta) + r1(2);
plot(xx, yy, 'r-', 'LineWidth', 1);
xlim([-3e-6 3e-6]); ylim([-3e-6 3e-6]);
theta = 0 : 0.01 : 2*pi;
radius = a;
xx = radius * cos(theta) + r2(1);
yy = radius * sin(theta) + r2(2);
plot(xx, yy, 'r-', 'LineWidth', 1);
xlim([-3e-6 3e-6]); ylim([-3e-6 3e-6]);
title('fy_per_vol');
colorbar;

fprintf('f1 = %g\n',f1);
fprintf('f2 = %g\n', f2);
end


fmag = f1(2);
            
%% Plot result
if(phi_trifold)
figure;
subplot(1,3,1);
pc = pcolor(XX/a,YY/a,phi); set(pc, 'EdgeColor', 'none');
title('\phi (Two Spheres)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,2);
pc = pcolor(XX/a,YY/a,phi_left); set(pc, 'EdgeColor', 'none');
title('\phi (Left Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,3);
pc = pcolor(XX/a,YY/a,phi_right); set(pc, 'EdgeColor', 'none');
title('\phi (Right Spheres)');
ylabel('y/a');
xlabel('x/a');
colorbar;
end



%% Hmag
if(hmag)
    absHXHY = sqrt(HXn.^2+HYn.^2);

figure;
pc = pcolor(XX/(2*a),YY/(2*a),absHXHY); 
set(pc, 'EdgeColor', 'none');
title('|H| (Two Spheres)');
ylabel('y/D');
xlabel('x/D');
colorbar; colormap('hot');
xlim([-2 2]);
ylim([-2 2]);

figure;
subplot(1,3,1);
pc = contourf(XX/(2*a),YY/(2*a),absHXHY, 20); 
% set(pc, 'EdgeColor', 'none');
title('|H| (Two Spheres)');
ylabel('y/D');
xlabel('x/D');
% rotate(pc, [0 0 1]', 90);
colorbar; colormap('hot');
% caxis([350 max(absHXHY,[],'all')]);
subplot(1,3,2);
pc = contourf(XX/(2*a),YY/(2*a),sqrt(HX_left.^2+HY_left.^2), 20); 
% set(pc, 'EdgeColor', 'none');
title('|H| (Left Sphere)');
ylabel('y/D');
xlabel('x/D');
% rotate(pc, [0 0 1]', 90);
colorbar; colormap('hot');
% caxis([350 max(sqrt(HX_left.^2+HY_left.^2),[],'all')]);
subplot(1,3,3);
pc = contourf(XX/(2*a),YY/(2*a),sqrt(HX_right.^2+HY_right.^2), 20); 
% set(pc, 'EdgeColor', 'none');
title('|H| (Right Sphere)');
ylabel('y/D');
xlabel('x/D');
% rotate(pc, [0 0 1]', 90);
colorbar; colormap('hot');
% caxis([350 max(sqrt(HX_right.^2+HY_right.^2),[],'all')]);
end

if(for_paper)
    absHXHY = sqrt(HXn.^2+HYn.^2);

    figure;
    pc = pcolor(XX/(2*a),YY/(2*a),absHXHY); 
    set(pc, 'EdgeColor', 'none');
    title('|H| (Two Grains)');
    ylabel('z/D');
    xlabel('s/D');
    colorbar; colormap('hot');
    axis equal;
    xlim([-2 2]);
    ylim([-2 2]);
    
    figure;
    pc = pcolor(XX/(2*a),YY/(2*a),perm_map_debug_two_sph/perm_free_space); 
    set(pc, 'EdgeColor', 'none');
    title('\mu/\mu_0 (Two Grains)');
    ylabel('z/D');
    xlabel('s/D');
    axis equal;
    xlim([-2 2]);
    ylim([-2 2]);
    colorbar;
end

%% HX Trifold
if(hx_trifold)
figure;
subplot(1,3,1);
pc = pcolor(XX/a,YY/a,HX/H0(1)); set(pc, 'EdgeColor', 'none');
title('H_x/H_{0x} (Two Spheres)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,2);
pc = pcolor(XX/a,YY/a,HX_left/H0(1)); set(pc, 'EdgeColor', 'none');
title('H_x/H_{0x} (Left Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,3);
pc = pcolor(XX/a,YY/a,HX_right/H0(1)); set(pc, 'EdgeColor', 'none');
title('H_x/H_{0x} (Right Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
end

if(hy_trifold)
figure;
subplot(1,3,1);
pc = pcolor(XX/a,YY/a,HYn/H0(2)); set(pc, 'EdgeColor', 'none');
title('H_y/H_{0y} (Two Spheres)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,2);
pc = pcolor(XX/a,YY/a,HY_left/H0(2)); set(pc, 'EdgeColor', 'none');
title('H_y/H_{0y} (Left Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,3);
pc = pcolor(XX/a,YY/a,HY_right/H0(2)); set(pc, 'EdgeColor', 'none');
title('H_y/H_{0y} (Right Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
end

if(hx_trifold_no_normalization)
figure;
subplot(1,4,1);
pc = pcolor(XX/(2*a),YY/(2*a),HXn); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('H_x (Two Spheres)');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([-0.5 0.5]);
colorbar; colormap('hot');
subplot(1,4,2);
pc = pcolor(XX/(2*a),YY/(2*a),HX_left); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('H_x (Left Sphere)');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([-0.5 0.5]);
colorbar; colormap('hot');
subplot(1,4,3);
pc = pcolor(XX/(2*a),YY/(2*a),HX_right); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('H_x (Right Sphere)');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([-0.5 0.5]);
colorbar; colormap('hot');
subplot(1,4,4);
pc = pcolor(XX/(2*a),YY/(2*a),HX_right+HX_left-HX); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('\Delta H_x');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([-0.5 0.5]);
colorbar; colormap('hot');
end

if(hy_trifold_no_normalization)
figure;
subplot(1,3,1);
pc = pcolor(XX/a,YY/a,HYn); set(pc, 'EdgeColor', 'none');
title('H_y (Two Spheres)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,2);
pc = pcolor(XX/a,YY/a,HY_left); set(pc, 'EdgeColor', 'none');
title('H_y (Left Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,3);
pc = pcolor(XX/a,YY/a,HY_right); set(pc, 'EdgeColor', 'none');
title('H_y (Right Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
end

if(perm_map_debug)
figure;
subplot(1,3,1)
pc = pcolor(XX/a,YY/a,perm_map_debug_two_sph/perm_free_space); 
set(pc, 'EdgeColor', 'none');
title('\mu/\mu_0 (Two Spheres)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,2)
pc = pcolor(XX/a,YY/a,perm_map_debug_left_sph/perm_free_space); 
set(pc, 'EdgeColor', 'none');
title('\mu/\mu_0 (Left Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
subplot(1,3,3)
pc = pcolor(XX/a,YY/a,perm_map_debug_right_sph/perm_free_space); 
set(pc, 'EdgeColor', 'none');
title('\mu/\mu_0 (Right Sphere)');
ylabel('y/a');
xlabel('x/a');
colorbar;
end

if(interp_match_du_paper)
%% Solve for differential in field per Du/Toffletto/Biswal paper
dHx = HX_right + HX_left - H0(1)*ones(size(HXn)) - HXn;
dHy = HY_left + HY_right - H0(2)*ones(size(HYn)) - HYn;
figure;
% [xxx, yyy] = meshgrid(linspace(-5*a, 5*a, n*10), linspace(-5*a, 5*a, m*10));
% dHy_interp = interp2(XX, YY, dHy, xxx, yyy, 'cubic');
% dHx_interp = interp2(XX, YY, dHx, xxx, yyy, 'cubic');
subplot(1,2,1);
pc = pcolor(XX./(2*a),YY./(2*a),dHy); set(pc, 'EdgeColor', 'none');
%caxis([-15 10]');
xlim([-4/2 4/2]);
ylim([-3/2 3/2]);
xlabel('x/D');
ylabel('y/D');
rotate(pc, [0 0 1]', 90);
title('\Delta H_x');
colorbar; colormap('hot');

subplot(1,2,2);
pc = pcolor(XX./(2*a),YY./(2*a),dHx); set(pc, 'EdgeColor', 'none');
%caxis([-10 10]');
title('\Delta H_y');
xlim([-4/2 4/2]);
ylim([-3/2 3/2]);
xlabel('x/D');
ylabel('y/D');
rotate(pc, [0 0 1]', 90);
colorbar; colormap('hot');
end

if(four_plots)
figure; subplot(2,2,1); 
pc = pcolor(XX,YY,sqrt(HXn.^2+HYn.^2)); set(pc, 'EdgeColor', 'none');
title('Two Spheres (ABS)');
colorbar;
subplot(2,2,2); pc = pcolor(XX,YY,sqrt(HX_left.^2+HY_left.^2)); 
set(pc, 'EdgeColor', 'none');
colorbar;
title('Left Sphere');
subplot(2,2,3); pc = pcolor(XX,YY,sqrt(HX_right.^2+HY_right.^2)); 
set(pc, 'EdgeColor', 'none');
colorbar;
title('Right Sphere');
subplot(2,2,4); pc = pcolor(XX,YY,ones(size(HXn))*sqrt(H0(1).^2+H0(2).^2)); 
set(pc, 'EdgeColor', 'none');
title('Far Field');

figure; subplot(2,2,1); 
pc = pcolor(XX,YY,HX); set(pc, 'EdgeColor', 'none');
title('Two Spheres (X)');
colorbar;
subplot(2,2,2); pc = pcolor(XX,YY,HX_left); set(pc, 'EdgeColor', 'none');
colorbar;
title('Left Sphere');
subplot(2,2,3); pc = pcolor(XX,YY,HX_right); set(pc, 'EdgeColor', 'none');
colorbar;
title('Right Sphere');
subplot(2,2,4); pc = pcolor(XX,YY,ones(size(HXn))*H0(1)); 
set(pc, 'EdgeColor', 'none');
title('Far Field');

figure; subplot(2,2,1); 
pc = pcolor(XX,YY,HY); set(pc, 'EdgeColor', 'none');
title('Two Spheres (Y)');
colorbar;
subplot(2,2,2); pc = pcolor(XX,YY,HY_left); set(pc, 'EdgeColor', 'none');
colorbar;
title('Left Sphere');
subplot(2,2,3); pc = pcolor(XX,YY,HY_right); set(pc, 'EdgeColor', 'none');
colorbar;
title('Right Sphere');
subplot(2,2,4); pc = pcolor(XX,YY,ones(size(HYn))*H0(2)); 
set(pc, 'EdgeColor', 'none');
title('Far Field');
end
