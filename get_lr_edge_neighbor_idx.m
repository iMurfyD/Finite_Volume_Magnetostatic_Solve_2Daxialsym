function nx = get_lr_edge_neighbor_idx(i,j,n,m)
    if(j == 1)
        nx = n*(i-1) + j +1;
    elseif(j==n)
        nx = n*(i-1) + j - 1;
    else
        error('wtf');
    end
end