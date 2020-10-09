function idx = indexf(i,j,k,syst)
idx = i+syst.l*(j-1)+(syst.l*syst.m)*(k-1);
end