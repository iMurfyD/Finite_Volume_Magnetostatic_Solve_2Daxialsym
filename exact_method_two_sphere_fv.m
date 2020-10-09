%% Exact sphere method, finite volume
% Ian DesJardin
% August 2020

 clear all; close all; clc;

%% Set up parameters
a = 1.4e-6; % Radius of sphere
sep = 2.2*a;
r1 = [0 0 0]';
r2 = [0 0 0]'; 
perm_free_space = 4*pi*1.00000000082e-7; % H*m^-1 
                                         % permiability of free space
                                         
susc = 0.96; % Suscepibility of material, arbitrary

perm = perm_free_space*(1+susc); % Linear media, eqn 6.30 Griffiths
H0 = [0.0 0.0 477]'; % A/m

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

%% Set up the grid (l x m x n 3D grid)
l = 100;
m = 100;
n = 500;
xdom = linspace(-10*a, 10*a, l);
ydom = linspace(-10*a, 10*a, m);
zdom = linspace(-15*a, 15*a, n);
dx = xdom(2)-xdom(1);
dy = ydom(2)-ydom(1);
dz = zdom(2)-zdom(1);
[XX,YY,ZZ] = meshgrid(xdom,ydom,zdom);

syst = struct('l',l,'m',m,'n',n,'a',a,'dx',dx,'dy',dy,'dz',dz,...
              'XX',XX,'YY',YY,'ZZ',ZZ,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,....
              'alpha', 10000);%0.2517);

%% Form FV Matrix
fprintf('Setting up Linear System\n');
[A,b, perm_map_debug_two_sph] = setup_system_sparse(syst);
% fprintf('Solving Linear System with \\ Operator\n');
% u = A\(b);
% fprintf('Done Solving Linear System with \\ Operator\n');


fprintf('Computing Initial Solution Guess for PCG\n');
% init_guess = repelem((-ydom'*H0(3)),l*m);
fprintf('Done Computing Initial Solution Guess for PCG\n');
fprintf('Computing Cholesky Preconditioners for PCG\n');
L = ichol(A);
fprintf('Done Computing Cholesky Preconditioners for PCG\n');
tol = 1e-10;
fprintf('Solving Linear System with PCG Method to %g tolerance\n',tol);
[u_numeric, flag,relres,iter,resvec] = pcg(A,b,tol,2000,L,L');
fprintf('Done Solving with PCG Method\n');
% phi = spread_1D_into_3D(u,syst);
phi_numeric = spread_1D_into_3D(u_numeric,syst);

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
for i = 1:length(xdom)
    for j = 1:length(ydom)
        for k = 1:length(zdom)
        idx = indexf(i,j,k,syst);
        rvec = [XX(i,j,k)-r1(1) YY(i,j,k)-r1(2) ZZ(i,j,k)-r1(3)]';
        r = norm(rvec,2);
        thta = atan2(rvec(1),rvec(3));
        if(r >= a) % Outside
        pt = (perm-perm_free_space)/(perm+2*perm_free_space);
        phi_theoretical(i,j) = -H0(3)*r*cos(thta) + ...
                           H0(3)*power(a,3)*pt*cos(thta)/power(r,2);
        else % Inside
            phi_theoretical(i,j) = -H0(3)*((3*perm_free_space)/...
                                           (perm+2*perm_free_space))*...
                                     r*cos(thta);
    
        end
        u_theoretical(idx) = phi_theoretical(i,j);
        end
    end
end

% [HX, HY, HZ] = gradient(phi,xdom,ydom,zdom);
% HX=-HX;
% HY=-HY;
[HXt, HYt, HZt] = gradient(phi_theoretical,xdom,ydom,zdom);
HXt=-HXt;
HYt=-HYt;

[HXn,HYn, HZn] = gradient(phi_numeric, xdom, ydom,zdom);
HXn=-HXn;
HYn=-HYn;

figure; pc=pcolor(squeeze(HZn(25,:,:)));set(pc,'EdgeColor','none');colorbar;


figure; subplot(2,3,1);
pc = pcolor(XX./a,YY./a,sqrt(HX.^2+HY.^2)); set(pc, 'EdgeColor', 'none');
colorbar; title('|H| Calculated');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,2);
pc = pcolor(XX./a,YY./a,sqrt(HXt.^2+HYt.^2)); set(pc, 'EdgeColor', 'none');
colorbar; title('|H| Theoretical');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,3);
pc = pcolor(XX./a, YY./a, abs(sqrt(HX.^2+HY.^2)-sqrt(HXt.^2+HYt.^2))./sqrt(HXt.^2+HYt.^2)); 
set(pc, 'EdgeColor', 'none');
colorbar; title('|\Delta |H| |/ |H|_t');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,4);
pc = pcolor(XX./a,YY./a,sqrt(HXn.^2+HYn.^2)); set(pc, 'EdgeColor', 'none');
colorbar; title('|H| Conjugate Gradient');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,5);
pc = pcolor(XX./a, YY./a, abs(sqrt(HXn.^2+HYn.^2)-sqrt(HXt.^2+HYt.^2))./sqrt(HXt.^2+HYt.^2)); 
set(pc, 'EdgeColor', 'none');
colorbar; title('|\Delta |H_n| |/ |H|_t');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,6);
pc = pcolor(XX./a, YY./a, abs(sqrt(HXn.^2+HYn.^2)-sqrt(HX.^2+HY.^2))./sqrt(HX.^2+HY.^2)); 
set(pc, 'EdgeColor', 'none');
colorbar; title('|\Delta |H_n - H| |/ |H|');
xlim([-2 2]); ylim([-2 2]);

figure; subplot(2,3,1);
pc = pcolor(XX./a,YY./a,HX); set(pc, 'EdgeColor', 'none');
colorbar; title('Hx Calculated');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,2);
pc = pcolor(XX./a,YY./a,HXt); 
set(pc, 'EdgeColor', 'none');
colorbar; title('Hx Theoretical');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,3);
pc = pcolor(XX./a, YY./a, (HX-HXt)); 
set(pc, 'EdgeColor', 'none');
colorbar; title('\Delta Hx');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,4);
pc = pcolor(XX./a,YY./a,HXn); set(pc, 'EdgeColor', 'none');
colorbar; title('Hxn Calculated');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,5);
pc = pcolor(XX./a,YY./a,HXn-HXt); 
set(pc, 'EdgeColor', 'none');
colorbar; title('Hxn-Hxt');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,6);
pc = pcolor(XX./a, YY./a, (HXn-HX)); 
set(pc, 'EdgeColor', 'none');
colorbar; title('Hxn-Hx');
xlim([-2 2]); ylim([-2 2]);

figure; subplot(2,3,1);
pc = pcolor(XX./a,YY./a,HY); set(pc, 'EdgeColor', 'none');
colorbar; title('Hy Calculated');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,2);
pc = pcolor(XX./a,YY./a,HYt);
set(pc, 'EdgeColor', 'none');
colorbar; title('Hy Theoretical');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,3);
pc = pcolor(XX./a, YY./a, (HY-HYt)./HYt); 
set(pc, 'EdgeColor', 'none');
colorbar; title('\Delta Hy / H_yt');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,4);
pc = pcolor(XX./a,YY./a,HYn); set(pc, 'EdgeColor', 'none');
colorbar; title('Hyn Calculated');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,5);
pc = pcolor(XX./a,YY./a,(HYn-HYt)./HYt); 
set(pc, 'EdgeColor', 'none');
colorbar; title('(Hyn-Hyt)/Hyt');
xlim([-2 2]); ylim([-2 2]);

subplot(2,3,6);
pc = pcolor(XX./a, YY./a, (HYn-HY)); 
set(pc, 'EdgeColor', 'none');
colorbar; title('Hyn-Hy');
xlim([-2 2]); ylim([-2 2]);

end    
[HX, HY] = gradient(phi,xdom,ydom);
HY=-HY;
HX=-HX;
% H = -grad(phi), HX is calculated with opposite sign convention I think

%% Form FV Matrix and solve for just left sphere
r2 = r1;
[A,b, perm_map_debug_left_sph] = setup_system_sparse(syst);
if(spy_mat)
figure;
spy([A b]);
title('Left Sphere A Matrix');
end
u = A\b;
% % u_guess = A\b;
% u = minres(A,b,tol, maxit, [], [], u_guess);
% % [u,flag,relres,iter,resvec] = gmres(A,b, [], tol, maxit, [], [], u_guess);
% fprintf('Finished second system\n');
% % disp(relres);
% % disp(resvec);
phi_left = spread_1D_into_2D(u, m,n);
[HX_left, HY_left] = gradient(phi_left,xdom,ydom);
HX_left = -HX_left;
HY_left = -HY_left;
% H = -grad(phi), HX is calculated with opposite sign convention I think

%% Form FV Matrix and solve for just right sphere
r2 = [0 sep/2]';
r1 = r2;
[A,b, perm_map_debug_right_sph] = setup_system_sparse(syst);
if(spy_mat)
figure;
spy([A b]);
title('Right Sphere A Matrix');
end
u = A\b;
% % u_guess = A\b;
% u = minres(A,b,tol, maxit, [], [], u_guess);
% % [u,flag,relres,iter,resvec] = gmres(A,b, [], tol, maxit, [], [], u_guess);
% fprintf('Finished third system\n');
% % disp(relres);
% % disp(resvec);
phi_right = spread_1D_into_2D(u, m,n);
[HX_right, HY_right] = gradient(phi_right,xdom,ydom);
HX_right=-HX_right;
HY_right=-HY_right;
% H = -grad(phi), HX is calculated with opposite sign convention I think


%% Formulate the maxwell stress tensor and integrate force around sphere
f1 = [0 0]';
f2 = [0 0]';
[FMX, FMY] = gradient(perm_map_debug_two_sph,xdom,ydom);
FMX = -0.5*(HX.^2+HY.^2).*FMX;
FMY = -0.5*(HX.^2+HY.^2).*FMY;

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
% fprintf('test_FV = %g \n', test_fV);
for ii = 1:m
    for jj = 1:n
        if(YY(ii,jj) <= 0) % r1
%             if(norm([XX(ii,jj) YY(ii,jj)]'-r1)<=a)
            if(abs(XX(ii,jj)) <= mean([dx dy]))
                circum = 2*pi*mean([dx dy])/2;
            else
                circum = 2*pi*abs(XX(ii,jj))/2;
            end
            c_map(ii,jj) = circum;
%             f_per_vol_dV = [0 test_fV]'*dx*dy*abs(circum);
            f_per_vol_dV = [FMX(ii,jj) FMY(ii,jj)]'*dx*dy*abs(circum);
            fx_per_vol_dV(ii,jj) = f_per_vol_dV(1);
            fy_per_vol_dV(ii,jj) = f_per_vol_dV(2);  
            f1 = f1 + f_per_vol_dV;
%             end
        else % r2
%             if(norm([XX(ii,jj) YY(ii,jj)]'-r2)<=a)
            if(abs(XX(ii,jj)) <= mean([dx dy]))
                circum = mean([dx dy]);
            else
                circum = 2*pi*abs(XX(ii,jj))/2;
            end
            c_map(ii,jj) = circum;
            f_per_vol_dV = [FMX(ii,jj) FMY(ii,jj)]'*dx*dy*abs(circum);
%             f_per_vol_dV = [0 test_fV]'*dx*dy*abs(circum);
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
    absHXHY = sqrt(HX.^2+HY.^2);

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
