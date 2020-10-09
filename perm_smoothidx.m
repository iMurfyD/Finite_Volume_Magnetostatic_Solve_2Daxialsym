function p = perm_smoothidx(idx,syst)
    [i,j,k] = indexfinv(idx,syst);
    r = [syst.XX(i,j,k), syst.YY(i,j,k), syst.ZZ(i,j,k)]';
    p = perm_smooth(r,syst);
end