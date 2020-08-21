function fz = calc_f_FV_parallel(sepz, Hymag, susc,a,n,m)
%calc_f_FV_parallel Calcualtes a bunch of cases at once
fz = zeros(size(sepz));
n = length(sepz);
parfor i = 1:n
    fz(i) = calc_truth_f(sepz(i), Hymag, susc, a,n,m);
end
end

