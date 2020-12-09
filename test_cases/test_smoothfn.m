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
n = 2000;
m = 2000;
xdom = linspace(-20*a, 20*a, n);
ydom = linspace(-20*a, 20*a, m);
ds = xdom(2)-xdom(1);
dz = ydom(2)-ydom(1);
[XX,YY] = meshgrid(xdom,ydom);

syst = struct('n',n,'m',m,'a',a,'ds',ds,'dz',dz,'XX',XX,'YY',YY,...
              'r1',r1,'r2',r2,'perm',perm,'pfs',perm_free_space,'H0',H0,....
              'alpha', 10000);%0.2517);
          
yz = ydom(floor(m/2+1):end);
smoofz = zeros(size(yz));

for idx = 1:length(yz)
    r = [0 yz(idx)]';
    smoofz(idx) = smoothfni(r,syst.r1,syst);
end

figure;
plot(yz/syst.a, smoofz, 'o-');
title('Smoothed Profile Along y-axis');
xlabel('y/a');
ylabel('\lambda(y)');
ylim([-0.1 1.1]);
xlim([0.5 1.5]);

