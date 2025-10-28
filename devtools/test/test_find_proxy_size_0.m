% Make sure things return with correct shape in 2D
addpath('../utils');
addpath('..');

rad = 10.0;
tol = 1e-13;
dim = 2;
half_sidelen = 0.5;

k = @(s,t) log_kernel(s,t);

[n_reg_pts, n_ring_pts] = find_proxy_size(@log_kernel, ...
    half_sidelen, dim, rad, tol);


dim = 3;

[n_reg_pts, n_ring_pts] = find_proxy_size(@log_kernel, ...
    half_sidelen, dim, rad, tol);