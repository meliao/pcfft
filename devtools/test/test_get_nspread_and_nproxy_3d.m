addpath(genpath('../../pcfft'));

% 3D
tol = 1e-8;
dim = 3;
halfside = 0.5;

k = @(s,t) one_over_r_kernel(s,t);

[grid_info, proxy_info] = get_nspread_and_nproxy(k, dim, tol, halfside);