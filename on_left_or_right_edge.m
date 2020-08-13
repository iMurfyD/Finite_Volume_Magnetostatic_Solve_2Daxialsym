function tf = on_left_or_right_edge(i,j,n,m)
    if(j == 1 || j ==n)
        tf = 1;
    else
        tf = 0;
    end
end