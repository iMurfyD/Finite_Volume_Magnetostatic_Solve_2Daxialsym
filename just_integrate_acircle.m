clear all; clc; close all;

% gridz = 500:1:1000;
gridz = 100:10:15000;
percent_err = zeros(size(gridz));

idx = 0;
for n = gridz

    idx = idx +1;
a = 1.4e-6;
da = 16;
xdom = linspace(-da*a, da*a, n);
dx = xdom(2)- xdom(1);
ydom = linspace(-da*a, da*a, n);
dy = ydom(2) - ydom(1);
[XX, YY] = meshgrid(xdom, ydom);
err = dx;

FF = zeros(size(XX));


for j = 1:n
    for i = 1:n
        r = abs(sqrt((XX(j,i)^2+YY(j,i)^2)));
     if((a-a/10) <= r && r <= a)
         FF(j,i) = 1;
     end
    end
end

circum = 0.0;

[TT,RR] = cart2pol(XX,YY);
FFmid = FF(round(n/2),:);
dt = (2*pi)/(sum(FF,'all')/(sum(FFmid)/2));
dr = (RR(2,2) - RR(2,3))*sqrt(2);%((1/sqrt(2))+1)/2;
% dr = a/max(diff(find(FFmid>0)));

for j = 1:n
    for i = 1:n
%         if FF(i,j) ~= 0
%             fprintf('neq\n');
%         end
% %         circ_val(i,j) = FF(i,j)*abs(RR(i,j))*dr*dt;
        circum = circum + FF(i,j)*abs(RR(i,j))*dr*dt;
    end
end

%     circum = trapz(dy.*trapz(FF, 2), 1);

truth = pi*(a^2 - (a-a/10)^2);
circum_val(idx) = circum;
percent_err(idx) = (circum_val(idx) - truth)/(truth);
% rat(idx) = circum_val(idx)/(2*pi*a);
if (mod(n, 1000) == 0)
    fprintf('Running grid = %g\n', n);
    fprintf('error = %g%% \n', abs(percent_err(idx))*100);
end

end

fido = fopen('circle_conv_data.txt', 'w');
for k1 = 1:length(gridz)
    fprintf(fido, '%.2f %.2f %.2f\n', gridz(k1), percent_err(k1), circum_val(k1));
end
fclose(fido);

% system('/usr/local/bin/gnuplot -p circle_integration_comparison.plt');

figure; 
plot(gridz, abs(percent_err)*100);%abs(percent_err)*100);
ylabel('Percent Error');
yyaxis right;
plot(gridz, circum_val);
% yline(8*pi);
ylabel('Integral');
xlabel('Grid Size (nxn)');
title('Error in Integrating [a/10 a] Disk Area of Circle, using Polar Method [-16a 16a], with a=1.4e-6');