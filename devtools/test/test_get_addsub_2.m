addpath(genpath("../../pcfft"));
close all;
clear;


rad = 2.0;


% Set up two source points
rng(4);
n_src = 400;
source_pts = rand(2,n_src) - 0.5;

% target points are close but not exactly = source points
n_targ = 317;
target_pts = rand(2,n_targ) - 0.5;
disp(size(target_pts));


src_weights = zeros(n_src,1);
src_weights(1) = 1.0;
src_weights = src_weights(:);

% Define the kernel
k = @(s,t) log_kernel(s,t);

K_src_to_target = log_kernel(struct('r',source_pts), struct('r',target_pts));

target_vals = K_src_to_target * src_weights;
n_nbr = 100; % 10000 points / 500 is approximately 20 boxes

src_info = struct;
src_info.r = source_pts;
targ_info = struct;
targ_info.r = target_pts;


tol = 1e-08;

[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the bins, sources, and targets.

% Plot the source points with blue dots
figure(1);
scatter(source_pts(1,:), source_pts(2,:), 100, 'b.');
hold on;
% Plot the target points with red dots
scatter(target_pts(1,:), target_pts(2,:), 100, 'r.');

% Plot the index of each bin at its center
for i = 0:grid_info.nbin(1) * grid_info.nbin(2) - 1
    bin_ctr = bin_center(i, grid_info);
    text(bin_ctr(1), bin_ctr(2), num2str(i), 'HorizontalAlignment', 'center');
end

% Draw a circle of radius 2 * proxy rad around center of box 0
box0_ctr = bin_center(0, grid_info);
ring = get_ring_points(100, 2 * proxy_info.radius, box0_ctr);
plot(ring(1,:), ring(2,:), 'k--');
% Plot the reg gridpoints with black x's
scatter(grid_info.r(1,:), grid_info.r(2,:), 100, 'kx');


[A_spread_s, sort_info_s] = get_spread(k, k, src_info, ...
    grid_info, proxy_info);

[A_spread_t, sort_info_t] = get_spread(k, k, targ_info, ...
    grid_info, proxy_info);

assert(all(~isnan(A_spread_s(:))));
assert(all(~isinf(A_spread_s(:))));
assert(all(~isnan(A_spread_t(:))));
assert(all(~isinf(A_spread_t(:))));

A_addsub = get_addsub(k, k, src_info, targ_info, grid_info, ...
    proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);

% % A_addsub = A_add - A_sub;

% K_grid2grid = log_kernel(grid_info, grid_info);
% % Put zeros on the diagonal of K_grid2grid
% for i = 1:size(K_grid2grid, 1)
%     K_grid2grid(i, i) = 0;
% end

% term1 = A_addsub * src_weights;
% disp("main: term1: ")
% disp(term1)

% % term2 = A_sub * src_weights;
% % disp("main: term2: ")
% % disp(term2)

% AKA = A_spread_t.' * K_grid2grid * A_spread_s;
% disp("main: AKA: ")
% disp( AKA)

% term3 = AKA * src_weights;
% disp("main: term3: ")
% disp(term3)

% evals_approx = term1 + term3;

% disp("main: evals_approx: ")
% disp(evals_approx)
% disp("main: target_vals: ")
% disp(target_vals)

% errors_at_target = max(abs(evals_approx - target_vals)) / max(abs(target_vals));
% disp("errors_at_target: " + num2str(errors_at_target));


% % Print out A_add, A_sub, and AKA for debugging
% disp("main: A_addsub: ")
% disp(full(A_addsub))
% % disp("main: A_sub: ")
% % disp(full(A_sub))
% % disp("main: AKA: ")
% % disp(full(AKA))


% assert(errors_at_target < tol);

