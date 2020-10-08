function p = perm_smoothidx(idx,syst)
    i = fix((idx-1)/syst.n)+1;
    j = rem(idx-1,syst.n)+1;
    r = [syst.XX(i,j), syst.YY(i,j)]';
    p = perm_smooth(r,syst);
end