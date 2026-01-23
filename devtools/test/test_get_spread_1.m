% Tests that get_spread returns a spreading matrix with the correct shape.
addpath(genpath("../../pcfft"));

rad = 2.0;
ntarg = 100;
target_pts = get_ring_points(ntarg, rad * sqrt(2));

% Set up a random source
rng(3);
n_src = 2;
half_side_len = 1.0;
source_pts = (rand(2,n_src) - 0.5) * half_side_len;
% Cube the uniform IID points so there is a non-uniform distribution
source_pts(1,:) = source_pts(1,:).^3;
source_pts(2,:) = source_pts(2,:).^3;
disp("source_pts:")
disp(source_pts)

src_weights = rand(n_src,1);
src_weights = src_weights(:);

src_info = struct;
src_info.r = source_pts;

% Make the "target points" the same as the source points.
targ_info = struct;
targ_info.r = source_pts;

K_src_to_target = log_kernel(struct('r',source_pts), struct('r',target_pts));
disp("K_src_to_target shape: ")
disp(size(K_src_to_target))
disp("src_weights shape:")
disp(size(src_weights))
target_vals = K_src_to_target * src_weights(:);

% Define the kernel
k = @(s,t) log_kernel(s,t);

tol = 1e-05;

n_nbr = 500; % 10000 points / 500 is approximately 20 boxes
[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);


disp("grid_info.ngrid: ")
disp(grid_info.ngrid)
disp("grid_info.dx: ")
disp(grid_info.dx)
disp("grid_info.Lbd:")
disp(grid_info.Lbd)
% disp(grid_info.r(:, 1:100))

% Plot the source points and the grid points
% scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
% hold on;
% scatter(src_info.r(1,:), src_info.r(2,:), 'ro');

A_spread = get_spread(k, k, src_info, grid_info, proxy_info);

disp("A_spread shape:")
disp(size(A_spread))

n_grid_pts = size(grid_info.r, 2);
n_src = size(src_info.r, 2);
expected_shape = [n_grid_pts, n_src];
disp("Expected shape of A_spread:")
disp(expected_shape)

% assert that A_spread has shape (n_grid_pts, n_src)
assert(all(isequal(size(A_spread), expected_shape)));
