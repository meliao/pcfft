% Chech that the output of dx_nproxy gives us the desired error levels.
addpath(genpath("../../pcfft"));
clear;
close all;


% Set up 100 target points on a ring of radius sqrt(2) * rad

% Set up a random source
rng(3);
n_src = 100;
side_len = 1.0;


% Define the kernel
k = @(s,t) one_over_r_kernel(s,t);

halfside = 0.5 * side_len;
crad = 2.5;

R = sqrt(3.0) * halfside;
rad = crad * R;

target_pts = get_sphere_points(100, sqrt(3.0) * crad * halfside);


% tol_vals = [1e-03];
tol_vals = logspace(-3, -8, 40);
n_tol_vals = size(tol_vals, 2);
error_vals = ones(n_tol_vals, 1);
n_reg_vals = ones(n_tol_vals, 1);
n_proxy_vals = ones(n_tol_vals, 1);
nbinpts_vals = ones(n_tol_vals, 1);
nspread_vals = ones(n_tol_vals, 1);
for i = 1:n_tol_vals
    tol = tol_vals(i);
    disp("Main: Working on tol " + num2str(tol));

    [grid_info, proxy_info] = dx_nproxy(k, 3, tol, halfside);
    
    nspread = grid_info.nspread;
    dx = grid_info.dx;

    xx = -halfside + dx / 2 + (0:nspread - 1) * dx;
    yy = xx;
    zz = xx;
    [X, Y, Z] = meshgrid(xx, yy, zz);
    X = permute(X,[3,1,2]);
    Y = permute(Y,[3,1,2]);
    Z = permute(Z,[3,1,2]);
    box_pts = [X(:).'; Y(:).'; Z(:).'];

    disp("Main: xx: " + num2str(xx(:)'))
    
    
    % Get source points sampled from the bin, which is 
    % [-nbinpts*dx/2, nbinpts*dx/2]^2
    nbinpts = grid_info.nbinpts;
    dx = grid_info.dx;
    bin_sidelen = nbinpts * dx;
    disp("Main: Final bin size: " + num2str(bin_sidelen));

    source_pts = (rand(3,n_src) - 0.5) * bin_sidelen;

    % Print out max and min x value of source_pts
    disp("Main: source_pts x min: " + num2str(min(source_pts(1, :))));
    disp("Main: source_pts x max: " + num2str(max(source_pts(1, :))));

    src_weights = rand(n_src,1);
    src_weights = src_weights(:);

    src_info = struct;
    src_info.r = source_pts;

    targ_info = struct;
    targ_info.r = source_pts;
    % disp(size(source_pts));
    % disp(size(src_weights));
    % assert(false);

    K_src_to_target = one_over_r_kernel(struct('r',source_pts), struct('r',target_pts));
    target_vals = K_src_to_target * src_weights(:);

    proxy_pts = proxy_info.r;
    K_src_to_proxy = one_over_r_kernel(struct('r',source_pts), struct('r',proxy_pts));
    proxy_vals = K_src_to_proxy * src_weights(:);
    
    % Solve the least squares problem
    rhs = proxy_vals;
    lhs = one_over_r_kernel(struct('r',box_pts), struct('r',proxy_pts));
    weights_reg = lsqminnorm(lhs, rhs, tol / 10);
    
    % Evaluate the approximation
    K_reg_to_target = one_over_r_kernel(struct('r',box_pts), struct('r',target_pts));
    target_vals_approx = K_reg_to_target * weights_reg(:);
    
    errors_at_target = max(abs(target_vals_approx(:) - target_vals(:))) / max(abs(target_vals));
    disp("Main: For tol " + num2str(tol) + ", observed error: " + num2str(errors_at_target));
    error_vals(i) = errors_at_target;
    n_proxy_vals(i) = proxy_info.n_points_total;
    nbinpts_vals(i) = grid_info.nbinpts;
    nspread_vals(i) = grid_info.nspread;

end
figure(1);
clf;
subplot(2,1,1);
plot(tol_vals(:), error_vals(:), '.-')
hold on
plot(tol_vals(:), tol_vals(:), '--')
xscale('log')
ylabel("Observed error")
yscale('log')
xlabel("Tolerance")
grid on;
subplot(2,1,2);
plot(tol_vals(:), n_proxy_vals(:), '.-');
hold on;
plot(tol_vals(:), nbinpts_vals(:), '.-');
plot(tol_vals(:), nspread_vals(:), '.-');
legend("nproxy", "nbinpts", "nspread");
grid on;
xscale('log');
yscale('log');

% figure(2);
% plot(target_vals_approx(:) - target_vals(:), '.-');
% % hold on;
% % plot(target_vals(:));
% % legend("Approx", "Exact");
% grid on;
% xlabel("Target point index");
% ylabel("Error in target eval");