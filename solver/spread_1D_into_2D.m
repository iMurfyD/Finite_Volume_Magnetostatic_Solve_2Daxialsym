function matu = spread_1D_into_2D(u,syst)
%spread_1D_into_2D Spread 1D array into 2D matrix following 
%convention used in rest of code
m = syst.m;
n = syst.n;
matu = zeros(m,n);
for i = 1:m
    for j = 1:n
        idx = n*(i-1) + j;
        matu(i,j) = u(idx);
    end
end
end

