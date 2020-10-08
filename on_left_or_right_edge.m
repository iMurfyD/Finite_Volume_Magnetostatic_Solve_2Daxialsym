function tf = on_left_or_right_edge(i,j,n,m)
%on_left_or_right_edge Returns 1 if on left or right edge of mesh, 0 otherwise
%
%     i = y coordinate on grid from bottom (i = 0 <=> y = 0)
%         1 <= i <= n  (MATLAB indexes from 1)
%     j = x coordinate on grid from left (j = 0 <=> x = 0)
%         1 <= i <= m  (MATLAB indexes from 1)
%         confusing since matricies are index (row column), not (x,y)
%     n = number of grid points in x direction (columns in matrix)
%     m = number of grid points in y direction (rows in matrix)

    if(j == 1 || j ==n)
        tf = 1;
    else
        tf = 0;
    end
end
