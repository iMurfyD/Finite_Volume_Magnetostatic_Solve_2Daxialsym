function x = smoothfni(r,ri,dx,dy, a)
%% Calculate magnetic permiabiliy with smoothing function
% Thickness ratio, tune so fn is continuous but 
% grains are still defined in a sphere like shape
thicc = 0.4*min(dx, dy); 
x= (1/2)*(tanh((a - norm(r-ri))/thicc)+1); 
end