function [A,b, perm_map_debug] = setup_system_sparse(m,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy,H0)
%Given grid information and location of spheres, creates
%linear system to solve for current free magnetic potential at each node
b = zeros(n*m, 1);
perm_map_debug = zeros(size(XX));
self_to_self = zeros(length(XX),1);
self_to_right = zeros(length(XX),1);
self_to_up = zeros(length(XX),1);
for i = 1:m
    for j = 1:n
        idx = n*(i-1) + j; % Ordering of nodes in 1D
        perm_map_debug(i,j) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy);
        if on_corner(i,j,n,m)
            [nx, ny] = get_corner_neighbor_idx(i,j,n,m);
            self_to_self(idx) = -(dy*perm_avg(idx,nx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                                -(dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
            if(i == 1 && j ==1) % Lower left, handle dot products accordingly
                self_to_right(idx) = dy*perm_avg(idx,nx, n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
                self_to_up(idx) = dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            elseif(i == 1 && j == n) % Lower right
                self_to_right(idx) = 0.0;
                self_to_up(idx) = dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            elseif(i == m && j == 1) % Top left
                self_to_up(idx) = 0.0;
                self_to_right(idx) = dy*perm_avg(idx,nx, n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            elseif(i == m && j == n) % Top right
                self_to_right(idx) = 0.0;
                self_to_up(idx) = 0.0;
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy + ...
                         perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            else
                error('Something went wrong');
            end
        elseif on_top_or_bot_edge(i,j,n,m)
            ny = get_topbot_neighbor_idx(i,j,n,m);
            self_to_right(idx) = dy*perm_avg(idx,idx+1,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            self_to_self(idx) = -(dy*perm_avg(idx,idx-1,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                                -(dy*perm_avg(idx,idx+1,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                                -(dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
            if(i == 1) % On bottom, dot(surface normal,H)/H = -1
                self_to_up(idx) = dx*perm_avg(idx,ny,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            else % On top, dot(surface normal,H) = 1 
                self_to_up(idx) = 0.0;
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(2)*dx;
            end
        elseif on_left_or_right_edge(i,j,n,m)       
            aboveme = n*(i-1) + j - n;
            belowme = n*(i-1) + j + n;
            nx = get_lr_edge_neighbor_idx(i,j,n,m);
            self_to_up(idx) = dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            self_to_self(idx) = -(dy*perm_avg(idx,nx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                                -(dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy) + ...
                                -(dx*perm_avg(idx,belowme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
            if(j == 1) % Left edge, dot(surface normal,H)/H = -1
                self_to_right(idx) = dy*perm_avg(idx,nx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
                b(idx) = -perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy;
            else % right edge, dot(surface normal, H)/H = 1
                self_to_right(idx) = 0.0;
                b(idx) = perm_smoothidx(idx,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)*H0(1)*dy;
            end
        else
            aboveme = n*(i-1) + j - n;
            belowme = n*(i-1) + j + n;
            left = n*(i-1) + j - 1;
            right = n*(i-1) + j + 1;
            self_to_up(idx) = dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy;
            self_to_right(idx) = dy*perm_avg(idx,right,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx;
            self_to_self(idx) = -(dy*perm_avg(idx,left,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                                -(dy*perm_avg(idx,right,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dx) + ...
                                -(dx*perm_avg(idx,belowme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy) + ...
                                -(dx*perm_avg(idx,aboveme,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)/dy);
        end
        
    end
end

A = spdiags([[self_to_up(n+1:end-n); self_to_up(1:n); zeros(n,1)] self_to_right self_to_self ...
             [0; self_to_right(1:end-1)] [zeros(n,1);self_to_up(n+1:end-n); self_to_up(1:n)]], ...
            [-n -1 0 1 n], m*n, m*n);

end

function p = perm_avg(i, j, n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy)
    t = perm_smoothidx(i,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy) + ...
        perm_smoothidx(j,n,XX,YY,r1,r2,a,perm,perm_free_space,dx,dy);
    p = t/2;
end

