function nx = get_lr_edge_neighbor_idx(i,j,n,m)
% get_lr_neighbor_idx Get index of neighbors in grid on left or right edge
%     Using indexing consistent with meshgrid and starting at the bottom
%     left and going up row-wise from left to right, bottom to top
%     Calculated for a given i and j (coordinates) and n and m (grid size)
%     what index the x and y neighbor is for a corner 
%
%     i = y coordinate on grid from bottom (i = 0 <=> y = 0)
%         1 <= i <= n  (MATLAB indexes from 1)
%     j = x coordinate on grid from left (j = 0 <=> x = 0)
%         1 <= i <= m  (MATLAB indexes from 1)
%         confusing since matricies are index (row column), not (x,y)
%     n = number of grid points in x direction (columns in matrix)
%     m = number of grid points in y direction (rows in matrix)

    if(j == 1)
        nx = n*(i-1) + j +1;
    elseif(j==n)
        nx = n*(i-1) + j - 1;
    else
        error('wtf');
    end
end
