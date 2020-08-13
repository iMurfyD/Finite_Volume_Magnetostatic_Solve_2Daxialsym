function tf = on_top_or_bot_edge(i,j,n,m)
    if(i == 1 || i == m)
        tf =  1;
    else
        tf =  0;
    end
end