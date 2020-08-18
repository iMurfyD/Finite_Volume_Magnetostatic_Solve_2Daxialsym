function [A,b, perm_map_debug] = setup_system(m,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy,H0)
%Given grid information and location of spheres, creates
%linear system to solve for current free magnetic potential at each node
A = zeros(n*m, n*m);
b = zeros(n*m, 1);
perm_map_debug = zeros(size(XX));
for i = 1:m
    for j = 1:n
        idx = n*(i-1) + j; % Ordering of nodes in 1D
        perm_map_debug(i,j) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy);
        if on_corner(i,j,n,m)
            [nx, ny] = get_corner_neighbor_idx(i,j,n,m);
            A(idx, nx) = dy*perm_avg(idx,nx, n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            A(idx, ny) = dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            A(idx, idx) = -(dy*perm_avg(idx,nx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                          -(dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
            if(i == 1 && j ==1) % Lower left, handle dot products accordingly
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            elseif(i == 1 && j == n) % Lower right
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            elseif(i == m && j == 1) % Top left
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            elseif(i == m && j == n) % Top right
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            else
                error('Something went wrong');
            end
        elseif on_top_or_bot_edge(i,j,n,m)
            ny = get_topbot_neighbor_idx(i,j,n,m);
            A(idx, ny) = dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            A(idx, idx-1) = dy*perm_avg(idx,idx-1,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            A(idx, idx+1) = dy*perm_avg(idx,idx+1,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            A(idx, idx) = -(dy*perm_avg(idx,idx-1,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                          -(dy*perm_avg(idx,idx+1,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                          -(dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
            if(i == 1) % On bottom, dot(surface normal,H)/H = -1
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            else % On top, dot(surface normal,H) = 1 
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            end
        elseif on_left_or_right_edge(i,j,n,m)
            nx = get_lr_edge_neighbor_idx(i,j,n,m);
            A(idx, nx) = dy*perm_avg(idx,nx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            aboveme = n*(i-1) + j - n;
            A(idx, aboveme) = dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            belowme = n*(i-1) + j + n;
            A(idx, belowme) = dx*perm_avg(idx,belowme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            A(idx, idx) = -(dy*perm_avg(idx,nx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                          -(dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy) + ...
                          -(dx*perm_avg(idx,belowme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
            if(j == 1) % Left edge, dot(surface normal,H)/H = -1
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy;
            else % right edge, dot(surface normal, H)/H = 1
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy;
            end
        else
            aboveme = n*(i-1) + j - n;
            belowme = n*(i-1) + j + n;
            left = n*(i-1) + j - 1;
            right = n*(i-1) + j + 1;
            A(idx, aboveme) = dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            A(idx, belowme) = dx*perm_avg(idx,belowme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            A(idx, left) = dy*perm_avg(idx,left,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            A(idx, right) = dy*perm_avg(idx,right,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            A(idx, idx) = -(dy*perm_avg(idx,left,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                          -(dy*perm_avg(idx,right,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                          -(dx*perm_avg(idx,belowme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy) + ...
                          -(dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
        end
        
    end
end
end

function p = perm_avg(i, j, n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)
    t = perm_smoothidx(i,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy) + ...
        perm_smoothidx(j,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy);
    p = t/2;
end

