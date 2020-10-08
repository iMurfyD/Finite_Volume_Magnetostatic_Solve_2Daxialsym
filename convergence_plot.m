
clear all; clc; close all;


nmz = [100 200 300 400 500 600 700 800 900 1000 1100 1250 1350 1500]';% 1750 2000 2500 3000 4000 5000]';

fv_nmz_sol = zeros(size(nmz));

sep = 2.1;
H0 = 1.0;
susc = 1.0;
a = 1.0;

for i = 1:length(nmz)
    fv_nmz_sol(i) = calc_truth_f_circum_method(sep, H0, susc, a, nmz(i), nmz(i));
end


figure;
plot(nmz, fv_nmz_sol, 'o-');
title('Convergence of Solution with nxn Grid');
xlabel('n');
ylabel('F'); 
