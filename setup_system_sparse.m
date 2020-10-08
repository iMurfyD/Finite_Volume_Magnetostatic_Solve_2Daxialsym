function [A,b, perm_map_debug] = setup_system_sparse(syst)
%Given grid information and location of spheres, creates
%linear system to solve for current free magnetic potential at each node
m = syst.m;
n  = syst.n;
perm_free_space = syst.pfs;
dx = syst.dx;
dy = syst.dy;
H0 = syst.H0;

% if(H0(1) ~= 0.0)
%     error('H0 is nonzero in x direction\n');
% end

b = zeros(n*m, 1);
perm_map_debug = zeros(size(syst.XX));
self_to_self = zeros(n*m, 1);
self_to_right = zeros(n*m, 1);
self_to_up = zeros(n*m, 1);
self_to_down = zeros(n*m, 1);
self_to_left = zeros(n*m, 1);
for i = 1:m
    for j = 1:n
        idx = n*(i-1) + j; % Ordering of nodes in 1D
        aboveme = idx + n;
        belowme = idx - n;
        left = idx - 1;
        right = idx + 1;
        perm_map_debug(i,j) = perm_smoothidx(idx,syst);
        if on_corner(i,j,n,m)
            if(i == 1 && j ==1) % Lower left, handle dot products accordingly
                self_to_right(idx) =  perm_at_interface(i,j,idx,right, syst);
                self_to_up(idx) =  perm_at_interface(i,j,idx,aboveme,syst);
                b(idx) = perm_free_space*H0(1)*dy+ ...
                         perm_at_top_bot_edge(idx,syst)*H0(2)*dx;
            elseif(i == 1 && j == n) % Lower right
                self_to_up(idx) =  perm_at_interface(i,j,idx,aboveme,syst);
                self_to_left(idx) =  perm_at_interface(i,j,idx,left,syst);
                b(idx) = -perm_free_space*H0(1)*dy + ...
                         perm_at_top_bot_edge(idx,syst)*H0(2)*dx;
            elseif(i == m && j == 1) % Top left
                self_to_down(idx) =  perm_at_interface(i,j,idx,belowme,syst);
                self_to_right(idx) =  perm_at_interface(i,j,idx,right,syst);
                b(idx) =  perm_free_space*H0(1)*dy + ...
                         -perm_at_top_bot_edge(idx,syst)*H0(2)*dx;
            elseif(i == m && j == n) % Top right
                self_to_left(idx) = perm_at_interface(i,j,idx,left,syst);
                self_to_down(idx) = perm_at_interface(i,j,idx,belowme,syst);
                b(idx) = -perm_free_space*H0(1)*dy+ ...
                         -perm_at_top_bot_edge(idx,syst)*H0(2)*dx;
            else
                error('Something went wrong');
            end
        elseif (i == 1 || i == m) % Top or bottom edge

            self_to_right(idx) = perm_at_interface(i,j,idx,right,syst);
            self_to_left(idx) = perm_at_interface(i,j,idx,left,syst);
            if(i == 1) % On bottom, dot(surface normal,H)/H = -1
                self_to_up(idx) =  perm_at_interface(i,j,idx,aboveme,syst);
                b(idx) = perm_at_top_bot_edge(idx,syst)*H0(2)*dx;
            else % On top, dot(surface normal,H) = 1 
                self_to_down(idx) =  perm_at_interface(i,j,idx,belowme,syst);
                b(idx) = -perm_at_top_bot_edge(idx,syst)*H0(2)*dx;
            end
        elseif (j == 1 || j == n) % Left or right edge
            self_to_down(idx) =  perm_at_interface(i,j,idx,belowme,syst);
            self_to_up(idx) =  perm_at_interface(i,j,idx,aboveme,syst);
            if(j == 1) % Left edge, dot(surface normal,H)/H = -1
                self_to_right(idx) = perm_at_interface(i,j,idx,right,syst);
                b(idx) = perm_free_space*H0(1)*dy;
            else % right edge, dot(surface normal, H)/H = 1
                self_to_left(idx) = perm_at_interface(i,j,idx,left,syst);
                b(idx) = -perm_free_space*H0(1)*dy;
            end
        else
            self_to_down(idx) = perm_at_interface(i,j,idx,belowme,syst);
            self_to_left(idx) = perm_at_interface(i,j,idx,left,syst);
            self_to_up(idx) = perm_at_interface(i,j,idx,aboveme,syst);
            self_to_right(idx) = perm_at_interface(i,j,idx,right,syst);
        end
        
        self_to_self(idx) = self_to_up(idx)+self_to_down(idx)+...
                            self_to_left(idx)+self_to_right(idx);
        
    end
end

znx = zeros(n,1); 

A = spdiags([[-self_to_down(n+1:end); znx] [-self_to_left(2:end); 0] self_to_self ...
             [0; -self_to_right(1:end-1)] [znx; -self_to_up(1:end-n)]], ...
            [-n -1 0 1 n], m*n, m*n);
        
end

function p = perm_at_top_bot_edge(idx,syst)
    i = fix((idx-1)/syst.n)+1;
    j = rem(idx-1,syst.n)+1;
    rc = [syst.XX(i,j) syst.YY(i,j)]';
    if(i == 1)
        r = rc - [0 syst.dy/2]';
        p =  perm_smooth(r,syst);
    elseif(i == syst.m)
        r = rc + [0 syst.dy/2]';
        p =  perm_smooth(r,syst);
    end
end

function p = perm_at_interface(i,j,idx1, idx2, syst)
    % Not really an average anymore - just finds mu at interface of FV
    % cells
    i1 = fix((idx1-1)/syst.n)+1;
    if(i1 ~= i)
        error('i wrong');
    end
    j1 = rem(idx1-1,syst.n)+1;
    if(j1 ~= j)
        error('j wrong');
    end
    x1 = [syst.XX(i1,j1), syst.YY(i1,j1)]';
    i2 = fix((idx2-1)/syst.n)+1;
    j2 = rem(idx2-1,syst.n)+1;
    x2 = [syst.XX(i2,j2), syst.YY(i2,j2)]';
    r_interface = (x1+x2)/2;
    p = perm_smooth(r_interface,syst);
end

