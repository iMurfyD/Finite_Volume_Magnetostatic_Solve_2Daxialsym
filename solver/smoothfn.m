
function x = smoothfn(r,syst)
if (syst.r1 == syst.r2) % If Spheres are overlapping
    x = smoothfni(r,syst.r1,syst);
else
    x = smoothfni(r,syst.r1,syst) + smoothfni(r,syst.r2,syst);
end
end