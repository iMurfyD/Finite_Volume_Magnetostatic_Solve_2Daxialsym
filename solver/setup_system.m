function [A,b,perm_map_debug] = setup_system_sparse(syst)
%Given grid information and location of spheres, creates
%linear system to solve for current free magnetic potential at each node
m = syst.m;
n  = syst.n;
pfs = syst.pfs;
ds = syst.ds;
dz = syst.dz;
H0 = syst.H0;

b = zeros(n*m, 1);
perm_map_debug = zeros(n*m,1);
self_to_self = zeros(n*m, 1);
self_to_right = zeros(n*m, 1);
self_to_up = zeros(n*m, 1);
self_to_down = zeros(n*m, 1);
self_to_left = zeros(n*m, 1);

parfor idx = 1:syst.n*syst.m % index
    [i, j] = indexfinv(idx, syst);
    x = syst.XX(j,i);
    
    % Find indicies of neighbor volume elements 
    above = indexf(i,j+1,syst);
    below = indexf(i,j-1,syst);
    left = indexf(i-1,j,syst);
    right = indexf(i+1,j,syst);
    perm_map_debug(idx) = perm_smoothidx(idx,syst);   
    R = abs(x);
    
    % Sets "R" the radius in a cylindrical coordinate system
    if(x>ds)
        Rright = R-ds/2;
        Rleft = R+ds/2;
    elseif(R <= ds) % Prevents singularity
        Rright = ds;
        Rleft = ds;
        R = ds;
    else
        Rright = R+ds/2;
        Rleft = R-ds/2;
    end

    
    % X direction fluxes
    if(i == 1) % Touching the left face
        self_to_right(idx) = perm_at_interface(idx,right,syst)*Rright*dz/ds;
        b(idx) = b(idx)-pfs*H0(1)*Rright*dz;
    elseif(i == n) % Touching right face
        self_to_left(idx) = perm_at_interface(idx,left,syst)*Rleft*dz/ds;
        b(idx) = b(idx)+pfs*H0(1)*Rleft*dz;
    else % Somewhere inbetween the boundaries
        self_to_left(idx) = perm_at_interface(idx,left,syst)*Rright*dz/ds;
        self_to_right(idx) = perm_at_interface(idx,right,syst)*Rleft*dz/ds;
    end
    
    % Y direction fluxes
    if(j == 1) % Touching the bottom face
        self_to_up(idx) = perm_at_interface(idx,above,syst)*R*ds/dz;
        b(idx) = b(idx)-pfs*H0(2)*R*ds;
    elseif(j == m) % Touching the top face
        self_to_down(idx) = perm_at_interface(idx,below,syst)*R*ds/dz;
        b(idx) = b(idx)+pfs*H0(2)*R*ds;
    else % Somewhere inbetween the boundaries
        self_to_up(idx) = perm_at_interface(idx,above,syst)*R*ds/dz;
        self_to_down(idx) = perm_at_interface(idx,below,syst)*R*ds/dz;
    end
    
    self_to_self(idx) = self_to_up(idx) + self_to_down(idx)+...
                        self_to_left(idx) + self_to_right(idx);

end

b = -b; % Needs to swap sign convention with self_to_self being postitive

ynx = zeros(n,1);


% Constructs tridiagonal sparse system matrix
% Each element in the vector passed to spdiages is a diagonal
% Done to save memory and make conputation feasible
A = spdiags([[-self_to_down(n+1:end); ynx] ...
             [-self_to_left(2:end); 0] ...
             self_to_self ...
             [0; -self_to_right(1:end-1)] ...
             [ynx; -self_to_up(1:end-n)]], ...
            [-n -1 0 1 n ], m*n, m*n);
        
% perm_map_debug = perm_map_debug';
        
end

function p = perm_at_interface(idx1, idx2, syst)
    % Not really an average anymore - just finds mu at interface of FV
    % cells
    [i1, j1] = indexfinv(idx1, syst);
    x1 = [syst.XX(j1,i1), syst.YY(j1,i1)]';    
    [i2, j2] = indexfinv(idx2, syst);
    x2 = [syst.XX(j2,i2), syst.YY(j2,i2)]';
    r_interface = (x1+x2)/2;
    p = perm_smooth(r_interface,syst);
end




