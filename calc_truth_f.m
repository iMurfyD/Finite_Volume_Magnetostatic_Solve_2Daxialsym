function fmag = calc_truth_f(sep, Hymag, susc,a,n,m)
%calc_truth_f Calculates force around sphere in two sphere problem
%    Solves for H field around sphere and integrates
%    div Maxwell stress tensor through volume of sphere
%
%    Method from:
%        Numerical calculation of interaction forces between paramagnetic colloids in two-dimensional systems
%        by Di Du, Frank Toffoletto, and Sibani Lisa Biswal (2014)
%        Physical Review E 89, 043306
%
%    sep = normalized separation (r/a)
%           where r = center to center distance between spheres
%                 a = radius of sphere
%           example: r/a = 2 <==> spheres touching
%    Hymag = magnitude of applied H field (A/m)
%    susc = magnetic susceptibility of spheres
%    a = radius of sphere (m)
%    n = number of cells used in y direction
%    m = number of cells used in x direction
%        yes this is backward thank you meshgrid


if(Hymag == 0.0) 
    error('Need to apply a field');
end

%% Set up parameters
% a = 1.4e-6; % Radius of sphere
sep = sep*a; % Unnormalize it
r1 = [0 -sep/2]';
r2 = [0 sep/2]';
perm_free_space = 4*pi*1.00000000082e-7; % H*m^-1 
                                         % permiability of free space
                                        
perm = perm_free_space*(1+susc); % Linear media, eqn 6.30 Griffiths
H0 = [0 Hymag]'; % A/m

%% Set up the grid (n x m 2D grid)
xdom = linspace(-3*sep, 3*sep, n);
ydom = linspace(-3*sep, 3*sep, m);
dx = xdom(2)-xdom(1);
dy = ydom(2)-ydom(1);
[XX,YY] = meshgrid(xdom,ydom);

%% Form FV Matrix
[A,b, perm] = setup_system_sparse(m,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy,H0);
u = A\b;
phi = spread_1D_into_2D(u, m,n);
[HX, HY] = gradient(phi,xdom,ydom);
% H = -grad(phi), HX is calculated with opposite sign convention I think
HX = HX;
HY = -HY;

%% Formulate the maxwell stress tensor and integrate force around sphere
f1 = [0 0]';
f2 = [0 0]';
[FMX, FMY] = gradient(perm.*perm./perm_free_space,xdom,ydom);
FMX = -0.5*(HX.^2+HY.^2).*FMX;
FMY = -0.5*(HX.^2+HY.^2).*FMY;

figure; 
pc = pcolor(XX,YY,FMX); set(pc, 'EdgeColor', 'none');
% xlim([-3 3]); ylim([-3 3]);
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
% xlim([-3 3]); ylim([-3 3]);
title('FMY');
colorbar;

% border_threshold = sqrt(dx^2+dy^2);
d_map = zeros(size(XX));
fx_per_vol_dV = zeros(size(XX));
fy_per_vol_dV = zeros(size(XX));
for ii = 1:m
    for jj = 1:n
        r = [XX(ii,jj) YY(ii,jj)]';
        if( YY(ii,jj)<0)%norm(r1-r) <= a)%+(dx+dy)/4)
            rel_vec = r - r1;
            if(norm(rel_vec) <= a+(mean([dx dy])/2))
                % Calculate depth at each corner of cell, average
                % Lower left
                x = rel_vec(1)-(dx/2);
                y = rel_vec(2)-(dy/2);
                if(norm([x y]) > a)
                    d1 = 0.0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d1 = 2*a*sin(theta)*sin(phi);
                end
                % Lower right
                x = rel_vec(1)+(dx/2);
                y = rel_vec(2)-(dy/2);
                if(norm([x y]) > a)
                    d2 = 0.0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d2 = 2*a*sin(theta)*sin(phi);
                end
                % Upper right
                x = rel_vec(1)+(dx/2);
                y = rel_vec(2)+(dy/2);
                if(norm([x y]) > a)
                    d3 = 0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d3 = 2*a*sin(theta)*sin(phi);
                end
                % Upper left
                x = rel_vec(1)-(dx/2);
                y = rel_vec(2)+(dy/2);
                if(norm([x y]) > a)
                    d4 = 0.0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d4 = 2*a*sin(theta)*sin(phi);
                end
                depth = (d1+d2+d3+d4)/4;
            else
                depth = 0;
            end
            d_map(ii,jj) = depth;
            f_per_vol_dV = [FMX(ii,jj) FMY(ii,jj)]'*dx*dy*abs(depth);
            fx_per_vol_dV(ii,jj) = f_per_vol_dV(1);
            fy_per_vol_dV(ii,jj) = f_per_vol_dV(2);            
            f1 = f1 + f_per_vol_dV;
        else%if( norm(r2-r) <=a)%+(dx+dy)/4)
            rel_vec = r - r2;
            if(norm(rel_vec) <= a+(mean([dx dy])/2))
                                % Calculate depth at each corner of cell, average
                % Lower left
                x = rel_vec(1)-(dx/2);
                y = rel_vec(2)-(dy/2);
                if(norm([x y]) > a)
                    d1 = 0.0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d1 = 2*a*sin(theta)*sin(phi);
                end
                % Lower right
                x = rel_vec(1)+(dx/2);
                y = rel_vec(2)-(dy/2);
                if(norm([x y]) > a)
                    d2 = 0.0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d2 = 2*a*sin(theta)*sin(phi);
                end
                % Upper right
                x = rel_vec(1)+(dx/2);
                y = rel_vec(2)+(dy/2);
                if(norm([x y]) > a)
                    d3 = 0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d3 = 2*a*sin(theta)*sin(phi);
                end
                % Upper left
                x = rel_vec(1)-(dx/2);
                y = rel_vec(2)+(dy/2);
                if(norm([x y]) > a)
                    d4 = 0.0;
                else
                    theta = acos(y/a);
                    phi = acos(x/(a*sin(theta)));
                    d4 = 2*a*sin(theta)*sin(phi);
                end
                depth = (d1+d2+d3+d4)/4;
            else
                depth = 0;
            end
            d_map(ii,jj) = depth;
            f_per_vol_dV = [FMX(ii,jj) FMY(ii,jj)]'*dx*dy*abs(depth);
            fx_per_vol_dV(ii,jj) = f_per_vol_dV(1);
            fy_per_vol_dV(ii,jj) = f_per_vol_dV(2);  
            f2 = f2 + f_per_vol_dV;
        end
    end
end
    
fmag = f1(2);

figure; 
pc = pcolor(XX,YY,d_map); set(pc, 'EdgeColor', 'none');
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
% xlim([-3 3]); ylim([-3 3]);
title('depth');
colorbar;   

figure; 
pc = pcolor(XX,YY,fx_per_vol_dV); set(pc, 'EdgeColor', 'none');
% xlim([-3 3]); ylim([-3 3]);
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
theta = 0 : 0.01 : 2*pi;
radius = a;
xx = radius * cos(theta) + r2(1);
yy = radius * sin(theta) + r2(2);
plot(xx, yy, 'r-', 'LineWidth', 1);
% xlim([-3 3]); ylim([-3 3]);
title('fy_per_vol');
colorbar;

disp(f1);
disp(f2);

end
