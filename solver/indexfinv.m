function [i1, j1] = indexfinv(idx1, syst)
j1 = idivide(idx1 - 1, int32(syst.n))+1;
i1 = idx1-1 - ((j1-1)*syst.n) + 1;
end