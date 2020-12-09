
clear all; clc; close all;
addpath('../solver/', '..');

nmz = 200:200:2000;
fv_nmz_sol = zeros(size(nmz));


sep = 2.2;
H0 = 477.0;
susc = 0.96;
a = 1.4e-6;


for i = 1:length(nmz)
    fv_nmz_sol(i,:) = calc_f_two_grain(sep, H0, susc, ...
                                       a, nmz(i), nmz(i));
end

figure;
plot(nmz, fv_nmz_sol, 'o-');
title('Convergence of \int f_y 16ax16a domain Polar integration');
xlabel('n');
ylabel('sum f_y');

