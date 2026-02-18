% Makes sure get_spread returns without error on a 2D input.
addpath(genpath('../../pcfft'));


tol = 1e-10;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 13;
rng(0);
src_info_2d.r = (rand(2, n_src) - 0.5);
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
ntarg = 17;
targ_info_2d.r = rand(2, ntarg) - 0.5;


[grid_info, proxy_info] = get_grid(k, ...
    src_info_2d, targ_info_2d, tol);

disp("grid_info.ngrid: ")
disp(grid_info.ngrid)
disp("grid_info.dx: ")
disp(grid_info.dx)
disp(grid_info.r(:, 1:100))

scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
hold on;
scatter(src_info_2d.r(1,:), src_info_2d.r(2,:), 'ro');

kern_0 = @(s,t) log_kernel(s,t);
kern = @(s,t) log_kernel(s,t);


A_spread = get_spread(kern_0, kern, src_info_2d, grid_info, proxy_info);

% Check the size of A_spread.
assert(all(size(A_spread) == [size(grid_info.r, 2) n_src]));
assert(all(~isnan(A_spread(:))));
assert(all(~isinf(A_spread(:))));

%% 3D case


tol = 1e-6;
dim = 3;

k = @(s,t) one_over_r_kernel(s,t);

src_info_3d = struct;
n_src = 13;
rng(0);
src_info_3d.r = (rand(3, n_src) - 0.5);
src_info_3d.weights = rand(n_src, 1);

targ_info_3d = struct;
ntarg = 17;
targ_info_3d.r = rand(3, ntarg) - 0.5;


[grid_info, proxy_info] = get_grid(k, ...
    src_info_3d, targ_info_3d, tol);

disp("grid_info.ngrid: ")
disp(grid_info.ngrid)
disp("grid_info.dx: ")
disp(grid_info.dx)
% disp(grid_info.r(:, 1:100))

% scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
% hold on;
% scatter(src_info_2d.r(1,:), src_info_2d.r(2,:), 'ro');

kern_0 = @(s,t) one_over_r_kernel(s,t);
kern = @(s,t) one_over_r_kernel(s,t);


A_spread = get_spread(kern_0, kern, src_info_3d, grid_info, proxy_info);

% Check the size of A_spread.
assert(all(size(A_spread) == [size(grid_info.r, 2) n_src]));
assert(all(~isnan(A_spread(:))));
assert(all(~isinf(A_spread(:))));