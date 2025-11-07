addpath(genpath("../../pcfft"));

rad = 2.0;
ntarg = 200;
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

src_weights = zeros(n_src,1);
src_weights(2) = 1.0;
src_weights = src_weights(:);

src_info = struct;
src_info.r = source_pts;

% Make the "target points" the same as the source points.
targ_info = struct;
targ_info.r = source_pts;

K_src_to_target = log_kernel(source_pts, target_pts);
disp("K_src_to_target shape: ")
disp(size(K_src_to_target))
disp("src_weights shape:")
disp(size(src_weights))
target_vals = K_src_to_target * src_weights;

% Define the kernel
k = @(s,t) log_kernel(s,t);

tol = 1e-07;

n_nbr = 500; % 10000 points / 500 is approximately 20 boxes
[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);


% Beef up proxy_info
% nproxy = proxy_info.n_points_total * 10;
% proxy_pts = get_ring_points(nproxy, proxy_info.radius);
% proxy_info.r = proxy_pts;
% proxy_info.n_points_total = nproxy;

nbin = 2;

disp("grid_info.ngrid: ")
disp(grid_info.ngrid)
disp("grid_info.dx: ")
disp(grid_info.dx)
disp("grid_info.Lbd:")
disp(grid_info.Lbd)

disp("grid_info.r shape:")
disp(size(grid_info.r))
% disp(grid_info.r(:, 1:100))

% Plot the source points and the grid points
figure(1);
scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
hold on;
scatter(src_info.r(1,:), src_info.r(2,:), 'ro');

A_spread = get_spread(k, k, src_info, grid_info, proxy_info, nbin);

disp("A_spread shape:")
disp(size(A_spread))

reg_weights = A_spread * src_weights;
reg_weights = full(reg_weights);

K_reg_to_target = log_kernel(grid_info.r, target_pts);
disp("K_reg_to_target shape: ")
disp(size(K_reg_to_target))
disp("reg_weights shape:")
disp(size(reg_weights))

approx_targ_vals = K_reg_to_target * reg_weights;

errors_at_target = max(abs(approx_targ_vals - target_vals));
disp("target_vals: " + num2str(target_vals'));
disp("approx_targ_vals: " + num2str(approx_targ_vals'));
disp("errors_at_target: " + num2str(errors_at_target));

figure(2);
plot(target_vals, 'bo-');
hold on;
plot(approx_targ_vals, 'rx-');
legend('Exact Target Values', 'Approx Target Values');
figure(3);
plot(target_vals - approx_targ_vals, 'k.-');
title('Errors at Target Points');