function fmag = calc_truth_f(sep, Hymag, susc,a,n,m)
%calc_truth_f Calculates force around sphere in two sphere problem
%   Solves for H field around sphere and integrates
%   Maxwell stress tensor

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
[FMX, FMY] = gradient(perm,xdom,ydom);
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
         if( YY(ii,jj)<0)%norm(r1-r) <= a)%+(dx+dy)/4)
            f1 = f1 + [FMX(ii,jj) FMY(ii,jj)]'*dx*dy;
         else%if( norm(r2-r) <=a)%+(dx+dy)/4)
            f2 = f2 + [FMX(ii,jj) FMY(ii,jj)]'*dx*dy;
         end
    end
end
    
fmag = f1(2);
         
end