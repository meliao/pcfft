% 2D_0
% Make sure things return with correct shape in 2D
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
n_targ = 17;
targ_info_2d.r = (rand(2, n_targ) - 0.5) * half_sidelen;



[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);

assert(grid_info.dim == dim);
assert(proxy_info.dim == dim);
assert(ndims(grid_info.r) == 2);
assert(size(grid_info.r, 1) == 2);
disp("test: ngrid: ")
disp(size(grid_info.ngrid))
disp("test: mods: ") 
disp(mod(grid_info.ngrid, grid_info.nbinpts))
assert(all(mod(grid_info.ngrid, grid_info.nbinpts) == 0))


%% 3D_0

dim = 3;
src_info_3d = struct;
n_src = 13;
src_info_3d.r = (rand(3, n_src) - 0.5) * half_sidelen;
src_info_3d.weights = rand(n_src, 1);

targ_info_3d = struct;
targ_info_3d.radius = 4.0;
targ_info_3d.r = (rand(3, n_targ) - 0.5) * half_sidelen;
tol = 1e-08;

[grid_info, proxy_info] = get_grid(@one_over_r_kernel, ...
    src_info_3d, targ_info_3d, tol);

assert(grid_info.dim == dim);
assert(proxy_info.dim == dim);


figure(1);clf
scatter3(src_info_3d.r(1,:), src_info_3d.r(2,:), src_info_3d.r(3,:))
hold on
scatter3(grid_info.r(1,:), grid_info.r(2,:), grid_info.r(3,:),'.')
hold off

assert(ndims(grid_info.r) == 2);
assert(size(grid_info.r, 1) == 3);
