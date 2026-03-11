addpath(genpath('../../pcfft'));

%% 2D test with log_kernel

dim = 2;
tol = 1e-6;
halfside = 0.5;

[dx, nspread, nbinpts, proxy_info] = dx_nproxy(@log_kernel, dim, tol, halfside);

% dx is consistent with nspread and halfside
assert(abs(dx - 2*halfside/nspread) < 1e-14);

% nbinpts is floor(nspread/2)
assert(nbinpts == floor(nspread/2));

% proxy_info fields are consistent with inputs
assert(proxy_info.dim == dim);
assert(proxy_info.tol == tol);
assert(proxy_info.halfside == halfside);

% proxy point array has correct shape
assert(size(proxy_info.r, 1) == dim);
assert(size(proxy_info.r, 2) == proxy_info.n_points_total);


%% 3D test with one_over_r_kernel

dim = 3;
tol = 1e-6;
halfside = 0.5;

[dx, nspread, nbinpts, proxy_info] = dx_nproxy(@one_over_r_kernel, dim, tol, halfside);

assert(abs(dx - 2*halfside/nspread) < 1e-14);
assert(nbinpts == floor(nspread/2));
assert(proxy_info.dim == dim);
assert(proxy_info.tol == tol);
assert(proxy_info.halfside == halfside);
assert(size(proxy_info.r, 1) == dim);
assert(size(proxy_info.r, 2) == proxy_info.n_points_total);


%% Tighter tolerance requires at least as many spreading points (2D)

dim = 2;
halfside = 0.5;
tol_coarse = 1e-4;
tol_fine   = 1e-8;

[~, nspread_coarse] = dx_nproxy(@log_kernel, dim, tol_coarse, halfside);
[~, nspread_fine]   = dx_nproxy(@log_kernel, dim, tol_fine,   halfside, 2, true);

assert(nspread_fine >= nspread_coarse);
