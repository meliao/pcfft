addpath(genpath('../../pcfft'));

% 3D
tol = 1e-8;
dim = 3;
halfside = 0.5;

k = @(s,t) one_over_r_kernel(s,t);

[dx, nspread, nbinpts, proxy_info] = dx_nproxy(k, dim, tol, halfside);