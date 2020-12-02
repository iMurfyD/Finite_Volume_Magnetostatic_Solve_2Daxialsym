function fmag = calc_truth_f_circum_method(sep, Hymag, susc,a,n,m)
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

debug = 0;
diagnostic = 1;

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

if(diagnostic)
fprintf('Allocating space for grid at %s\n', datestr(now,'HH:MM:SS.FFF'));
end

%% Set up the grid (n x m 2D grid)
sdom = linspace(-16*a, 16*a, n);
zdom = linspace(-16*a, 16*a, m);
ds = sdom(2)-sdom(1);
dz = zdom(2)-zdom(1);
[XX,YY] = meshgrid(sdom,zdom);

syst = struct('m',m,'n',n,'a',a,'ds',ds,'dz',dz,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,....
              'alpha', 0.2517);
if(diagnostic)
fprintf('Finished allocating space for grid at %s\n', datestr(now,'HH:MM:SS.FFF'));
end

%% Form FV Matrix
if(diagnostic)
fprintf('Creating connectivity matrix at %s\n', datestr(now,'HH:MM:SS.FFF'));
end
[A,b, permmdbg] = setup_system_sparse(syst);
if(diagnostic)
fprintf('Finished creating connectivity matrix at %s\n', datestr(now,'HH:MM:SS.FFF'));
fprintf('Solving for scalar potential at %s\n', datestr(now,'HH:MM:SS.FFF'));
end
u = A\b;
phi = spread_1D_into_2D(u, syst);
perm = spread_1D_into_2D(permmdbg, syst);
if(diagnostic)
fprintf('Finished scalar potential at %s\n', datestr(now,'HH:MM:SS.FFF'));
fprintf('Solving for field at %s\n', datestr(now,'HH:MM:SS.FFF'));
end
[HX, HY] = gradient(phi,ds,dz);
HX(abs(HX) < 10*eps) = 0.0;
HY(abs(HY) < 10*eps) = 0.0;
HX = -HX;
HY = -HY;
if(diagnostic)
fprintf('Finished field at %s\n', datestr(now,'HH:MM:SS.FFF'));
end

if(debug)
    figure;
    pc = pcolor(XX./a,YY./a,HX); set(pc, 'EdgeColor', 'none');
    xlim([-3 3]); ylim([-3 3]);
    title('HX');
    xlabel('X');
    ylabel('Y');
    colorbar;
    
    figure;
    pc = pcolor(XX./a,YY./a,HY); set(pc, 'EdgeColor', 'none');
    xlim([-3 3]); ylim([-3 3]);
    title('HY');
    xlabel('X');
    ylabel('Y');
    colorbar;
    
    figure;
    pc = pcolor(XX./a,YY./a,sqrt(HX.^2+HY.^2)); set(pc, 'EdgeColor', 'none');
    xlim([-3 3]); ylim([-3 3]);
    title('|H|');
    xlabel('X');
    ylabel('Y');
    colorbar;
end
    
%% Formulate the maxwell stress tensor and integrate force around sphere
if(diagnostic)
fprintf('Creating Maxwell stress tensor at %s\n', datestr(now,'HH:MM:SS.FFF'));
end
f1 = [0 0]';
f2 = [0 0]';
[FMX, FMY] = gradient(perm,ds,dz);
if(debug)    
figure; 
pc = pcolor(XX./a,YY./a,perm); set(pc, 'EdgeColor', 'none');
xlim([-3 3]); ylim([-3 3]);
title('\mu');
colorbar; 

figure; 
pc = pcolor(XX./a,YY./a,FMX); set(pc, 'EdgeColor', 'none');
xlim([-3 3]); ylim([-3 3]);
title('\nabla \mu_x', 'interpreter', 'latex');
colorbar; 

figure; 
pc = pcolor(XX./a,YY./a,FMY); set(pc, 'EdgeColor', 'none');
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
xlim([-3 3]); ylim([-3 3]);
title('\nabla \mu_y');
colorbar;
end

FMX = -0.5*(HX.^2+HY.^2).*FMX;
FMY = -0.5*(HX.^2+HY.^2).*FMY;
if(diagnostic)
fprintf('Finished creating Maxwell stress tensor at %s\n', datestr(now,'HH:MM:SS.FFF'));
end
if(debug)    
figure; 
pc = pcolor(XX./a,YY./a,FMX); set(pc, 'EdgeColor', 'none');
xlim([-3 3]); ylim([-3 3]);
title('FMX');
colorbar; 

figure; 
pc = pcolor(XX./a,YY./a,FMY); set(pc, 'EdgeColor', 'none');
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
xlim([-3 3]); ylim([-3 3]);
title('FMY');
colorbar;
end

if(debug)
    fprintf('sep = %g\n', sep);
    fprintf('a = %g\n', a);
    fprintf('H0 = %g\n', Hymag);
    fprintf('susc = %g\n', susc);
    fprintf('n = %g\n', length(sdom));
    fprintf('m = %g\n', length(zdom));
end

% border_threshold = sqrt(dx^2+dy^2);
if(diagnostic)
fprintf('Integrating force at %s\n', datestr(now,'HH:MM:SS.FFF'));
end
[~,topidx] = min(abs(YY(:,1)-(r1(2)+a+3*dz)));
[~,botidx] = min(abs(YY(:,1)-(r1(2)-a-3*dz)));



circum = pi*abs(XX());
fy_per_vol_dV = FMY.*circum;

% Convert to polar coordinates centered around the test sphere
% Only uses subdomain ([-16a 16a] X, [botidx topidx] Y)
[TT,RR] = cart2pol(XX(botidx:topidx,:)-r1(1),YY(botidx:topidx,:)-r1(2));
fy_per_vol_dV_sd = fy_per_vol_dV(botidx:topidx,:); % slice of FMY in subdomain of interest
num_nonzeros_total = nnz(fy_per_vol_dV_sd);
num_nonzeros_midrow = 2;%nnz(fy_per_vol_dV_sd(round((topidx-botidx)/4),:));
dr = (RR(2,2) - RR(2,3))*sqrt(2);
[nzi, nzj, nzf] = find(fy_per_vol_dV_sd);
nztheta = zeros(num_nonzeros_total,1);
for k=1:num_nonzeros_total
    nztheta(k) = TT(nzi(k), nzj(k));
end

% Define angular bins to sort nonzero points into
num_theta_bins = num_nonzeros_total/2;
theta_bins = linspace(-pi, pi, num_theta_bins);

% Sort nonzero points into angle bins
nzbin = discretize(nztheta, theta_bins);

% find number of nonzeros in each angular bin (actually set to 1 if there
% are points)
num_points_in_bin = zeros(num_theta_bins,1);
for k = 1:num_nonzeros_total
    num_points_in_bin(nzbin(k)) = 1;%num_points_in_bin(nzbin(k)) + 1;
end


% For bins with zero members, assign theta-age to next nonzero
% bins above and below it
for b = 1:num_theta_bins
    if(num_points_in_bin(b) == 0)
        if(b == 1)
            % search down for next nonzero bin
            bdown = num_theta_bins;
            while(num_points_in_bin(bdown) == 0)
                bdown = bdown - 1;
            end
            % give that bin, this bins theta
            num_points_in_bin(bdown) = num_points_in_bin(bdown)/1.5;
            
            % search up for next nonzero bin
            bup = 2;
            while(num_points_in_bin(bup) == 0)
                bup = bup + 1;
            end
            % give that bin, this bins theta
            num_points_in_bin(bup) = num_points_in_bin(bup)/1.5;
        elseif(b == num_theta_bins)
            % search down for next nonzero bin
            bdown = num_theta_bins-1;
            while(num_points_in_bin(bdown) == 0)
                bdown = bdown - 1;
            end
            % give that bin, this bins theta
            num_points_in_bin(bdown) = num_points_in_bin(bdown)/1.5;
            
            % search up for next nonzero bin
            bup = 1;
            while(num_points_in_bin(bup) == 0)
                bup = bup + 1;
            end
            % give that bin, this bins theta
            num_points_in_bin(bup) = num_points_in_bin(bup)/1.5;
        else % one of the middle bins
            % search down for next nonzero bin
            bdown = b-1;
            while(num_points_in_bin(bdown) == 0)
                if(bdown == 1)
                    bdown = num_theta_bins;
                else
                    bdown = bdown - 1;
                end
            end
            % give that bin, this bins theta
            num_points_in_bin(bdown) = num_points_in_bin(bdown)/1.5;
            
            % search up for next nonzero bin
            bup = b+1;
            while(num_points_in_bin(bup) == 0)
                if(bup == num_theta_bins)
                    bup = 1;
                else
                    bup = bup + 1;
                end
            end
            % give that bin, this bins theta
            num_points_in_bin(bup) = num_points_in_bin(bup)/1.5;
        end
    end
end

% also could add another step here which would check how many
% layers are in each bin and create an effective number of cell 
% thickness

% Assign how much dt points in a given bin get
dt_per_bin = (2*pi/num_theta_bins)^(-1) .*num_points_in_bin;
dt_per_bin = 1./dt_per_bin;

nzdt = zeros(num_nonzeros_total,1);
for k = 1:num_nonzeros_total
    nzdt(k) = dt_per_bin(nzbin(k));
end
    


fmag = 0;
for k = 1:num_nonzeros_total
    fmag = fmag + nzf(k)*abs(RR(nzi(k),nzj(k)))*dr*nzdt(k);
end

if(diagnostic)
fprintf('Finsihed integrating force at %s\n', datestr(now,'HH:MM:SS.FFF'));
end

if(debug)
fprintf('mean c = %g\n', mean(circum,'all'));
fprintf('mean c/(2*pi) = %g\n', mean(circum,'all')/(2*pi));


figure; 
pc = pcolor(XX./a,YY./a,c_map); set(pc, 'EdgeColor', 'none');
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
xlim([-3 3]); ylim([-3 3]);
title('depth');
colorbar; 
caxis([0 pi*a]);

figure; 
pc = pcolor(XX./a,YY./a,fx_per_vol_dV); set(pc, 'EdgeColor', 'none');
xlim([-3 3]); ylim([-3 3]);
title('fx_per_vol');
colorbar; 

figure; 
pc = pcolor(XX./a,YY./a,fy_per_vol_dV); set(pc, 'EdgeColor', 'none');
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
xlim([-3 3]); ylim([-3 3]);
title('fy_per_vol');
colorbar;

fprintf('f1 = %g\n',f1);
fprintf('f2 = %g\n', f2);
end

end
