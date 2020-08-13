function tf = on_corner(i,j,n,m)
    if(i==1 && j==1 || ...
       i==1 && j==n || ...
       i==m && j==1 || ...
       i==m && j==n)
        tf = 1;
    else
        tf = 0;
    end
end