function p = perm_smooth(r,r1,r2,a,perm,pfs,dx,dy)
        s = smoothfn(r,r1,r2,dx,dy,a);
        p = 1/((s/perm)+((1-s)/pfs));
end