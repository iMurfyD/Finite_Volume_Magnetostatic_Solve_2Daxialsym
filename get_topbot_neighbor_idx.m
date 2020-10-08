function ny = get_topbot_neighbor_idx(i,j,n,m)
% get_topbot_neighbor_idx Get index of neighbors in grid on top or bottom edge
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

    if(i == 1) 
        ny = n*(i-1) + j + n;
    elseif(i == m)
        ny = n*(i-1) + j - n;
    else
        error('wtf');
    end
end
