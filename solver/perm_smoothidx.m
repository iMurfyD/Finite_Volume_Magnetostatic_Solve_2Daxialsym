function p = perm_smoothidx(idx,syst)
    [i,j] = indexfinv(idx,syst);
    r = [syst.XX(j,i), syst.YY(j,i)]';
    p = perm_smooth(r,syst);
end