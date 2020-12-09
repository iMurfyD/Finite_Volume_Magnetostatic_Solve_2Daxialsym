function fz = calc_f_FV_parallel(sepz, Hymag, susc,a,n,m)
%calc_f_FV_parallel Calcualtes a bunch of cases at once
%    Uses the parfor command in matlab to iterate through 
%    the finite volume force calculation for a given set of 
%    input parameters. Relies on MATLAB local profile being
%    setup to use multiples cores on your machine. For an
%    educational copy of MATLAB 2020a, I had to change this
%    manually
%
%    sepz = vector of normalized separations (r/a)
%           where r = center to center distance between spheres
%                 a = radius of sphere
%           example: r/a = 2 <==> spheres touching
%    Hymag = magnitude of applied H field (A/m)
%    susc = magnetic susceptibility of spheres
%    a = radius of sphere (m)
%    n = number of cells used in y direction
%    m = number of cells used in x direction
%        yes this is backward - thank you meshgrid

fz = zeros(size(sepz));
num = length(sepz);
parfor i = 1:num
    fz(i) = calc_f_two_grain(sepz(i), Hymag, susc, a,n,m);
end
end

