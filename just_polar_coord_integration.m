r1 = [0 0];
a = 1;
n = 1000;
xdom = linspace(-10*a, 10*a, n);
ydom = linspace(-10*a, 10*a, n);


rdom = linspace(0, r1(2)+a, 500);
dr = rdom(2)-rdom(1);
tdom = linspace(0, 2*pi, 500);
dt = tdom(2)-tdom(1);
[THETA, R] = meshgrid(tdom ,rdom);
[Xpol, Ypol] = pol2cart (THETA, R);
FMYpol = interp2(XX,YY,FMY,Xpol,Ypol);

circum = pi*abs(XX());
circumpol = interp2(XX,YY,circum,Xpol,Ypol);
fy_per_vol_dV = FMYpol.*circumpol;

fmag = 2*pi*sum(sum(dt*dr.*R.*fy_per_vol_dV));