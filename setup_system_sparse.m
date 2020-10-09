function [A,b, perm_map_debug] = setup_system_sparse(syst)
%Given grid information and location of spheres, creates
%linear system to solve for current free magnetic potential at each node
m = syst.m;
n  = syst.n;
l = syst.l;
pfs = syst.pfs;
dx = syst.dx;
dy = syst.dy;
dz = syst.dz;
H0 = syst.H0;

b = zeros(n*m*l, 1);
perm_map_debug = zeros(size(syst.XX));
self_to_self = zeros(n*m*l, 1);
self_to_right = zeros(n*m*l, 1);
self_to_up = zeros(n*m*l, 1);
self_to_down = zeros(n*m*l, 1);
self_to_left = zeros(n*m*l, 1);
self_to_front = zeros(n*m*l, 1);
self_to_behind = zeros(n*m*l, 1);

for i = 1:l
for j = 1:m
for k = 1:n
    idx = indexf(i,j,k,syst); % Ordering of nodes in 1D
    above = indexf(i,j,k+1,syst);
    below = indexf(i,j,k-1,syst);
    left = indexf(i,j-1,k,syst);
    right = indexf(i,j+1,k,syst);
    front = indexf(i+1,j,k,syst);
    behind = indexf(i-1,j,k,syst);
    perm_map_debug(i,j) = perm_smoothidx(idx,syst);

    % X direction fluxes
    if(i == 1) % Touching the back face
        self_to_front(idx) = perm_at_interface(idx,front,syst)*dy*dz/dx;
        b(idx) = -pfs*H0(1)*dy*dz;
    elseif(i == l) % Touching front face
        self_to_behind(idx) = perm_at_interface(idx,behind,syst)*dy*dz/dx;
        b(idx) =  pfs*H0(1)*dy*dz;
    else % Somewhere inbetween the boundaries
        self_to_front(idx) = perm_at_interface(idx,front,syst)*dy*dz/dx;
        self_to_behind(idx) = perm_at_interface(idx,behind,syst)*dy*dz/dx;
    end
    
    % Y direction fluxes
    if(j == 1) % Touching the left face
        self_to_right(idx) = perm_at_interface(idx,right,syst)*dx*dz/dy;
        b(idx) = -pfs*H0(2)*dx*dz;
    elseif(j == m) % Touching the right face
        self_to_left(idx) = perm_at_interface(idx,left,syst)*dx*dz/dy;
        b(idx) =  pfs*H0(2)*dx*dz;
    else % Somewhere inbetween the boundaries
        self_to_right(idx) = perm_at_interface(idx,right,syst)*dx*dz/dy;
        self_to_left(idx) = perm_at_interface(idx,left,syst)*dx*dz/dy;
    end
    
    % Z direction fluxes
    if(k == 1) % Touching the bottom face
        self_to_up(idx) = perm_at_interface(idx,above,syst)*dx*dy/dz;
        b(idx) = -pfs*H0(3)*dx*dy;
    elseif(k == n) % Touching the top face
        self_to_down(idx) = perm_at_interface(idx,below,syst)*dx*dy/dz;
        b(idx) =  pfs*H0(3)*dx*dy;
    else
        self_to_up(idx) = perm_at_interface(idx,above,syst)*dx*dy/dz;
        self_to_down(idx) = perm_at_interface(idx,below,syst)*dx*dy/dz;
    end        

    self_to_self(idx) = self_to_up(idx)+self_to_down(idx)+...
                        self_to_left(idx)+self_to_right(idx) + ...
                        self_to_front(idx)+self_to_behind(idx);

end
end
end

% b = -b; % Needs to swap sign convention with self_to_self being postitive

znx = zeros(l*m,1);
ynx = zeros(l,1);


A = spdiags([[-self_to_down((l*m)+1:end); znx] ...
             [-self_to_left(l+1:end); ynx] ...
             [-self_to_behind(2:end); 0] ...
             self_to_self ...
             [0; -self_to_front(1:end-1)] ...
             [ynx; -self_to_right(1:end-l)]...
             [znx; -self_to_up(1:end-(l*m))]], ...
            [-l*m -l -1 0 1 l l*m], m*n*l, m*n*l);
        
end

function p = perm_at_interface(idx1, idx2, syst)
    % Not really an average anymore - just finds mu at interface of FV
    % cells
    [i1, j1, k1] = indexfinv(idx1, syst);
    x1 = [syst.XX(i1,j1,k1), syst.YY(i1,j1,k1), syst.ZZ(i1,j1,k1)]';    
    [i2, j2, k2] = indexfinv(idx2, syst);
    x2 = [syst.XX(i2,j2,k2), syst.YY(i2,j2,k2), syst.ZZ(i2,j2,k2)]';
    r_interface = (x1+x2)/2;
    p = perm_smooth(r_interface,syst);
end




