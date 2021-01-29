function x = smoothfni(r,ri,syst)
%% Calculate magnetic permiabiliy with smoothing function
% Thickness ratio, tune so fn is continuous but 
% grains are still defined in a sphere like shape
% Smoothed function approximation of discontinuous magnetic susceptibility
% Outlined in PHYSICAL REVIEW E 89, 043306 (2014)

thicc = syst.alpha*min(syst.ds, syst.dz); 
r = r/syst.a;
ri = ri/syst.a;
tanharg = (1 - norm(r-ri,2))/(thicc);
% tanharg = tanharg + 1;
x= (1/2)*(tanh(tanharg)+1);
end