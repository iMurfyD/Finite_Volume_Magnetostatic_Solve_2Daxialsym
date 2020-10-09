function [i1, j1, k1] = indexfinv(idx1, syst)
k1 = idivide(idx1-1,int32(syst.l*syst.m))+1;
j1 = idivide(idx1 - 1 - ((k1-1)*syst.l*syst.m), int32(syst.l))+1;
i1 = idx1-1 - ((k1-1)*syst.l*syst.m)- ((j1-1)*syst.l) + 1;
end