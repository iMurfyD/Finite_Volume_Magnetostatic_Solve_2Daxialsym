
function x = smoothfn(r,r1,r2,dx,dy,a)
if (r1 == r2) % If Spheres are overlapping
    x = smoothfni(r, r1,dx,dy,a);
else
    x = smoothfni(r, r1,dx,dy,a) + smoothfni(r, r2, dx, dy,a);
end
end