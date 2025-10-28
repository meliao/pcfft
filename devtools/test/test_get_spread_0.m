% Makes sure get_spread returns without error on a 2D input.
addpath(genpath('../../pcfft'));


rad = 10.0;
tol = 1e-13;
dim = 2;
half_sidelen = 0.5;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 13;
src_info_2d.r = (rand(2, n_src) - 0.5) * half_sidelen;
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
targ_info_2d.radius = 4.0;


[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);


kern_0 = @(s,t) log_kernel(s,t);
kern = @(s,t) log_kernel(s,t);


[A_spread, bin_info] = get_spread(kern_0, kern, src_info_2d, grid_info, proxy_info);

% Check the size of A_spread.
assert(all(size(A_spread) == [grid_info.ngrid n_src]));
assert(all(~isnan(A_spread(:))));
assert(all(~isinf(A_spread(:))));