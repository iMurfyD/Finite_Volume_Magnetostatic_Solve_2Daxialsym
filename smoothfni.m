function x = smoothfni(r,ri,syst)
%% Calculate magnetic permiabiliy with smoothing function
% Thickness ratio, tune so fn is continuous but 
% grains are still defined in a sphere like shape
thicc = syst.alpha*syst.a;%min(syst.dx, syst.dy); 
r = r/syst.a;
ri = ri/syst.a;
tanharg = (1 - norm(r-ri,2))/(thicc);
% tanharg = tanharg + 1;
x= (1/2)*(tanh(tanharg)+1);
end