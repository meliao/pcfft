addpath(genpath('../../pcfft/'));


rad = 10.0;
tol = 1e-12;
dim = 2;
halfside = 0.5;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 13;
src_info_2d.r = (rand(2, n_src) - 0.5) * halfside;
src_info_2d.weights = rand(n_src, 1);

[grid_info, proxy_info] = compute_nspread_nproxy(k, dim, tol, halfside);

