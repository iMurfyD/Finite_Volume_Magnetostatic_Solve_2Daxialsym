function [nx, ny] = get_corner_neighbor_idx(i,j,n,m)
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