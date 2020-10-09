function matu = spread_1D_into_3D(u,syst)
%spread_1D_into_2D Spread 1D array into 3D matrix following 
%convention used in rest of code
l = syst.l;
m = syst.m;
n = syst.n;
matu = zeros(l,m,n);
for i = 1:l
    for j = 1:m
        for k = 1:n
            idx = indexf(i,j,k,syst);
            matu(i,j,k) = u(idx);
        end
    end
end

