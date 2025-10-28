addpath('../utils');
% ring of radius <rad> discretized with <n_proxy_pts>
rad = 3.0;

target_pts = sqrt(2.0) * rad * [1.0 1.0 1.0].';

% Set up a random source
rng(3);
n_src = 100;
half_side_len = 1.0;

source_pts = (rand(3,n_src) - 0.5) * half_side_len;

src_weights = rand(n_src,1);
src_weights = src_weights(:);
% disp(size(source_pts));
% disp(size(src_weights));
% assert(false);

K_src_to_target = one_over_r_kernel(source_pts, target_pts);
target_vals = K_src_to_target * src_weights(:);

% Define the kernel
k = @(s,t) one_over_r_kernel(s,t);

tol_vals = [1e-04 1e-05 1e-06 1e-07 1e-08];
n_tol_vals = size(tol_vals, 2);
error_vals = ones(n_tol_vals, 1);
n_reg_vals = ones(n_tol_vals, 1);
for i = 1:n_tol_vals
    tol = tol_vals(i);
    % Get the number of discretization points necessary
    [n_reg_pts, n_proxy_pts] = find_proxy_size(k, half_side_len, 3, ...
        rad, tol);
    
    reg_pts = get_regular_grid(n_reg_pts , half_side_len, 3);
    proxy_pts = get_cube_points(n_proxy_pts , rad);
    K_src_to_proxy = one_over_r_kernel(source_pts, proxy_pts);
    % disp(size(K_src_to_proxy))
    proxy_vals = K_src_to_proxy * src_weights(:);
    
    % Solve the least squares problem
    rhs = proxy_vals;
    lhs = one_over_r_kernel(reg_pts, proxy_pts);
    weights_reg = lhs \ rhs;
    
    % Evaluate the approximation
    K_reg_to_target = one_over_r_kernel(reg_pts, target_pts);
    target_vals_approx = K_reg_to_target * weights_reg(:);
    
    errors_at_target = abs(target_vals_approx - target_vals);
    error_vals(i) = errors_at_target;
    n_reg_vals(i) = n_reg_pts;

end
disp(n_reg_vals)
plot(tol_vals(:), error_vals(:), '.-')
hold on
plot(tol_vals(:), tol_vals(:), '--')
xscale('log')
ylabel("Observed error")
yscale('log')
xlabel("Tolerance")

