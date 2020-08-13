function p = perm_smoothidx(idx, n,XX,YY,r1,r2,a,perm,pfs,dx,dy)
    i = fix((idx-1)/n)+1;
    j = rem(idx-1,n)+1;
    r = [XX(i,j), YY(i,j)]';
    p = perm_smooth(r,r1,r2,a,perm,pfs,dx,dy);
end