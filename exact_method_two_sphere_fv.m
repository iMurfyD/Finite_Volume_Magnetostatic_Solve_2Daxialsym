%% Exact sphere method, finite volume
% Ian DesJardin
% August 2020

% clear all; close all; clc;

%% Set up parameters
a = 1; % Radius of sphere
sep = 6*a;
r1 = [0 -sep/2]';
r2 = [0 sep/2]'; 
perm_free_space = 4*pi*1.00000000082e-7; % H*m^-1 
                                         % permiability of free space
                                         
susc = 1.0; % Suscepibility of material, arbitrary

perm = perm_free_space*(1+susc); % Linear media, eqn 6.30 Griffiths
H0 = [0 1]'; % A/m

%% Define what stuff to plot
four_plots = 0;
interp_match_du_paper = 1   ;
perm_map_debug=0;
spy_mat = 0;
hy_trifold = 0;
hx_trifold = 0;
hx_trifold_no_normalization =0;
hy_trifold_no_normalization=0;
phi_trifold=0;
hmag=1;

%% Set up the grid (n x m 2D grid)
n = 1000;
m = 1000;
xdom = linspace(-3*sep, 3*sep, n);
ydom = linspace(-3*sep, 3*sep, m);
dx = xdom(2)-xdom(1);
dy = ydom(2)-ydom(1);
[XX,YY] = meshgrid(xdom,ydom);

%% Form FV Matrix
[A,b, perm_map_debug_two_sph] = setup_system_sparse(m,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy,H0);
if(spy_mat)
figure;
spy([A b]);
title('Two Sphere A Matrix');
end
u = A\b;
phi = spread_1D_into_2D(u, m,n);
[HX, HY] = gradient(phi,xdom,ydom);
% H = -grad(phi), HX is calculated with opposite sign convention I think
HX = HX;
HY = -HY;

%% Form FV Matrix and solve for just left sphere
r2 = r1;
[A,b, perm_map_debug_left_sph] = setup_system_sparse(m,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy,H0);

if(spy_mat)
figure;
spy([A b]);
title('Left Sphere A Matrix');
end
u = A\b;
phi_left = spread_1D_into_2D(u, m,n);
[HX_left, HY_left] = gradient(phi_left,xdom,ydom);
% H = -grad(phi), HX is calculated with opposite sign convention I think
HX_left = HX_left;
HY_left = -HY_left;

%% Form FV Matrix and solve for just right sphere
r2 = [0 sep/2]';
r1 = r2;
[A,b, perm_map_debug_right_sph] = setup_system_sparse(m,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy,H0);

if(spy_mat)
figure;
spy([A b]);
title('Right Sphere A Matrix');
end
u = A\b;
phi_right = spread_1D_into_2D(u, m,n);
[HX_right, HY_right] = gradient(phi_right,xdom,ydom);
% H = -grad(phi), HX is calculated with opposite sign convention I think
HX_right = HX_right;
HY_right = -HY_right;



%% Formulate the maxwell stress tensor and integrate force around sphere
f1 = [0 0]';
f2 = [0 0]';
[FMX, FMY] = gradient(perm_map_debug_two_sph,xdom,ydom);
FMX = -0.5*(HX.^2+HY.^2).*FMX;
FMY = -0.5*(HX.^2+HY.^2).*FMY;


% border_threshold = sqrt(dx^2+dy^2);
for ii = 1:m
    for jj = 1:n
         r = [XX(ii,jj) YY(ii,jj)]';
%         if( norm(r2-r) <=a+border_threshold && ...
%             norm(r2-r) >= a-border_threshold) % Edge of sphere
%             nunit = (r2-r)/norm(r2-r);
%             max_stress_tens = zeros(2,2);
%             Hx = HX(ii,jj);
%             Hy = HY(ii,jj);
%             Hmag = sqrt(Hx^2+Hy^2);
%             max_stress_tens(1,1) = perm*Hx^2 - (1/2)*Hmag^2;
%             max_stress_tens(1,2) = perm*Hx*Hy;
%             max_stress_tens(2,1) = perm*Hy*Hx;
%             max_stress_tens(2,2) = perm*Hy^2 - (1/2)*Hmag^2;
%             f = f + max_stress_tens*nunit*dx;
%         end
        if( norm(r1-r) <=a+(dx+dy)/4 )
            f1 = f1 + [FMX(ii,jj) FMY(ii,jj)]'*dx*dy;
        end
        if( norm(r2-r) <=a+(dx+dy)/4 )
            f2 = f2 + [FMX(ii,jj) FMY(ii,jj)]'*dx*dy;
        end
    end
end
            
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
figure;
subplot(1,3,1);
absHXHY = sqrt(HX.^2+HY.^2);
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
pc = pcolor(XX/a,YY/a,HY/H0(2)); set(pc, 'EdgeColor', 'none');
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
pc = pcolor(XX/(2*a),YY/(2*a),HX); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('H_x (Two Spheres)');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([0.0 0.5]);
colorbar; colormap('hot');
subplot(1,4,2);
pc = pcolor(XX/(2*a),YY/(2*a),HX_left); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('H_x (Left Sphere)');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([0.0 0.5]);
colorbar; colormap('hot');
subplot(1,4,3);
pc = pcolor(XX/(2*a),YY/(2*a),HX_right); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('H_x (Right Sphere)');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([0.0 0.5]);
colorbar; colormap('hot');
subplot(1,4,4);
pc = pcolor(XX/(2*a),YY/(2*a),HX_right+HX_left-HX); set(pc, 'EdgeColor', 'none');
rotate(pc, [0 0 1]', 90);
title('\Delta H_x');
ylabel('y/D');
xlabel('x/D');
xlim([-2/4 0]);
ylim([0.0 0.5]);
colorbar; colormap('hot');
end

if(hy_trifold_no_normalization)
figure;
subplot(1,3,1);
pc = pcolor(XX/a,YY/a,HY); set(pc, 'EdgeColor', 'none');
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
dHx = HX_right + HX_left - H0(1)*ones(size(HX)) - HX;
dHy = HY_left + HY_right - H0(2)*ones(size(HY)) - HY;
figure;
[xxx, yyy] = meshgrid(linspace(-5*a, 5*a, n*10), linspace(-5*a, 5*a, m*10));
dHy_interp = interp2(XX, YY, dHy, xxx, yyy, 'cubic');
dHx_interp = interp2(XX, YY, dHx, xxx, yyy, 'cubic');
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
pc = pcolor(XX,YY,sqrt(HX.^2+HY.^2)); set(pc, 'EdgeColor', 'none');
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
subplot(2,2,4); pc = pcolor(XX,YY,ones(size(HX))*sqrt(H0(1).^2+H0(2).^2)); 
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
subplot(2,2,4); pc = pcolor(XX,YY,ones(size(HX))*H0(1)); 
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
subplot(2,2,4); pc = pcolor(XX,YY,ones(size(HY))*H0(2)); 
set(pc, 'EdgeColor', 'none');
title('Far Field');
end
