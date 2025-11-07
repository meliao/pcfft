addpath(genpath("../../pcfft"));

rad = 2.0;
target_pts = sqrt(2.0) * rad * [1.0 1.0].';

% Set up a random source
rng(3);
n_src = 1000;
half_side_len = 1.0;

source_pts = (rand(2,n_src) - 0.5) * half_side_len;

src_weights = rand(n_src,1);
src_weights = src_weights(:);

src_info = struct;
src_info.r = source_pts;

targ_info = struct;
targ_info.r = source_pts;
% disp(size(source_pts));
% disp(size(src_weights));
% assert(false);

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
    % Get the number of discretization points necessary. Set n_nbr = n_src
    % so that there is only one spreading box
    [grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_src);
    
    reg_pts = get_regular_grid(grid_info.ngrid(1), half_side_len, 2);
    proxy_pts = get_ring_points(proxy_info.n_points_total , rad);
    K_src_to_proxy = log_kernel(source_pts, proxy_pts);
    proxy_vals = K_src_to_proxy * src_weights(:);
    
    % Solve the least squares problem
    rhs = proxy_vals;
    lhs = log_kernel(reg_pts, proxy_pts);
    weights_reg = lhs \ rhs;
    
    % Evaluate the approximation
    K_reg_to_target = log_kernel(reg_pts, target_pts);
    target_vals_approx = K_reg_to_target * weights_reg(:);
    
    errors_at_target = abs(target_vals_approx - target_vals);
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

