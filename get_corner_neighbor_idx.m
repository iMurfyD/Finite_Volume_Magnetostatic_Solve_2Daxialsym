function [nx, ny] = get_corner_neighbor_idx(i,j,n,m)
% get_corner_neighbor_idx Get index of neighbors in grid from corner nodes
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

    if(i == 1 && j ==1) % Lower left
        nx = 2; ny = n+1;
    elseif(i == 1 && j == n) % Lower right
        nx = n-1; ny = 2*n;
    elseif(i == m && j == 1) % Top left
        ny = n*(m-2)+1; nx = (n*m) - (n-2);
    elseif(i == m && j == n) % Top right
        ny = (n*m)-n; nx = (n*m)-1;
    else
        error('What');
    end
end
