addpath(genpath("../../pcfft"));

target_rad = 4.01;
ntarg = 200;
target_pts = get_sphere_points(ntarg, target_rad);



% Set up two source points
rng(4);
n_src = 200;
rng(0);
source_pts = (rand(3, n_src) - 0.5) * 2;


src_weights = randn(n_src,1);
src_weights = src_weights(:);

% Define the kernel
k = @(s,t) one_over_r_kernel(s,t);

K_src_to_target = k(struct('r',source_pts), struct('r',target_pts));

target_vals = K_src_to_target * src_weights;
n_nbr = 500; % 10000 points / 500 is approximately 20 boxes



%% Part 2: Re-do the above but loop over tol values and plot tol vs error

tol_vals = [1e-02 1e-03 1e-04 1e-05 1e-06 1e-07 1e-08];
n_tol_vals = size(tol_vals, 2);
error_vals = zeros(n_tol_vals, 1);
dx_vals = zeros(n_tol_vals, 1);
nbinpts_vals = zeros(n_tol_vals, 1);
nproxy_vals = zeros(n_tol_vals, 1);

for i = 1:n_tol_vals
    tol = tol_vals(i);
    src_info = struct;
    src_info.r = source_pts;

    % Make the "target points" the same as the source points.
    targ_info = struct;
    targ_info.r = source_pts;
    [grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);


    % Loop through the bins and make sure that the bin center is > proxy_info.radius
    % away from the target_pts.
    n_bins_total = grid_info.nbin(1) * grid_info.nbin(2) * grid_info.nbin(3) ;
    for j = 0:n_bins_total -1
        [pts, center, row_idxes] = grid_pts_for_box_2d(j, grid_info);
        xdists = center(1) - target_pts(1);
        ydists = center(2) - target_pts(2);
        dists = sqrt(xdists.^2 + ydists.^2);
        assert(all(dists > proxy_info.radius));
    end

    A_spread = get_spread(k, k, src_info, grid_info, proxy_info);
    reg_weights = A_spread * src_weights;
    reg_weights = full(reg_weights);

    K_reg_to_target = one_over_r_kernel(grid_info, struct('r',target_pts));
    % disp("K_reg_to_target shape: ")
    % disp(size(K_reg_to_target))
    % disp("reg_weights shape:")
    % disp(size(reg_weights))

    approx_targ_vals = K_reg_to_target * reg_weights;

    errors_at_target = max(abs(approx_targ_vals - target_vals)) / max(abs(target_vals));

    % Save the error and dx vals
    error_vals(i) = errors_at_target;
    dx_vals(i) = grid_info.dx;
    nbinpts_vals(i) = grid_info.nbinpts;
    nproxy_vals(i) = proxy_info.n_points_total;
end

%%
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
plot(tol_vals(:), nproxy_vals(:), '.-');
hold on;
plot(tol_vals(:), nbinpts_vals(:), '.-');
% plot(tol_vals(:), nspread_vals(:), '.-');
legend("nproxy", "nbinpts");
grid on;
xscale('log');
