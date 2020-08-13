function ny = get_topbot_neighbor_idx(i,j,n,m)
    if(i == 1) 
        ny = n*(i-1) + j + n;
    elseif(i == m)
        ny = n*(i-1) + j - n;
    else
        error('wtf');
    end
end