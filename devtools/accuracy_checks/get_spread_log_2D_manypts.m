addpath(genpath("../../pcfft"));

rad = 2.0;
ntarg = 100;
target_pts = get_ring_points(ntarg, rad * sqrt(2));

% Set up a random source
rng(3);
n_src = 10000;
half_side_len = 1.0;
source_pts = (rand(2,n_src) - 0.5) * half_side_len;
% Cube the uniform IID points so there is a non-uniform distribution
source_pts(1,:) = source_pts(1,:).^3;
source_pts(2,:) = source_pts(2,:).^3;

src_weights = rand(n_src,1);
src_weights = src_weights(:);

src_info = struct;
src_info.r = source_pts;

% Make the "target points" the same as the source points.
targ_info = struct;
targ_info.r = source_pts;

K_src_to_target = log_kernel(source_pts, target_pts);
target_vals = K_src_to_target * src_weights(:);

% Define the kernel
k = @(s,t) log_kernel(s,t);

tol_vals = [1e-04 1e-05 1e-06 1e-07 1e-08];
n_tol_vals = size(tol_vals, 2);
error_vals = ones(n_tol_vals, 1);
n_reg_vals = ones(n_tol_vals, 1);
for i = 1:n_tol_vals
    tol = tol_vals(i);

    n_nbr = 500; % 10000 points / 500 is approximately 20 boxes
    [grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);

    nbin = grid_info.nspread;

    A_spread = get_spread(k, k, src_info, grid_info, proxy_info, nbin);

    reg_weights = A_spread * src_weights;
    reg_weights = full(reg_weights);

    K_reg_to_target = log_kernel(grid_info.r, target_pts);

    approx_targ_vals = K_reg_to_target * reg_weights;
    
    errors_at_target = max(abs(target_vals_approx - target_vals));
    error_vals(i) = errors_at_target;
    n_reg_vals(i) = grid_info.ngrid(1);

end
disp(n_reg_vals)
plot(tol_vals(:), error_vals(:), '.-')
hold on
plot(tol_vals(:), tol_vals(:), '--')
xscale('log')
ylabel("Observed error")
yscale('log')
xlabel("Tolerance")

