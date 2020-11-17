
clear all; clc; close all;

nmz = 100:200:2000;
% nmz = [100 200 300 400 500 600 700 800 900 1000]';
% nmz = [1500 1750 2000 2250 2500 2575 3000 3500 4000 4500 5000]';

% fv_nmz_sol = zeros(size(nmz));
fv_top_nmz_sol = zeros(size(nmz));


sep = 2.2;
H0 = 477.0;
susc = 0.96;
a = 1.4e-6;


for i = 1:length(nmz)
    fv_nmz_sol(i,:) = calc_truth_f_circum_method(sep, H0, susc, ...
                                               a, nmz(i), nmz(i));
%     fv_top_nmz_sol(i) = calc_f_top_circum_method(sep, H0, susc, ...
%                                                a, nmz(i), nmz(i));
end

figure;
plot(nmz, fv_nmz_sol, 'o-');
% plot(nmz, fv_nmz_sol, 'o-', nmz,sum(fv_nmz_sol,2));%, nmz);%, -fv_top_nmz_sol, '.');
title('Convergence of \int f_y 4ax4a domain Simpsons rule');
% title('Convergence of Solution with nxn Grid');
xlabel('n');
% legend('Bottom', 'Top', 'Top+Bot');
% ylabel('F');
ylabel('sum f_y');

