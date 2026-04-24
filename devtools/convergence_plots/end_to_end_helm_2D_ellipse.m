addpath(genpath("../../pcfft"));
close all;
clear;


% Set up an ellipse of sources, which will also be the targets.
rng(1);
n_src = 500;

src_pts = get_ring_points(n_src, 1.0);
% Shift by [0.3, 0.3]
src_pts(1,:) = src_pts(1,:) + 0.3;
src_pts(2,:) = src_pts(2,:) + 0.3;
% Scale it so it's aspect ratio 2:1
src_pts(1,:) = 2*src_pts(1,:);

% Target points are random in [-0.5, 0.5] x [-0.5, 0.5]
rng(2);
target_pts = rand(2, n_src) - 0.5;

% Plot the source points
figure(1);
clf;
scatter(src_pts(1,:), src_pts(2,:), 'x');
hold on;
scatter(target_pts(1,:), target_pts(2,:), 'o');
xlabel("x");
ylabel("y");

%% Compute distances between source points
dists = sqrt((src_pts(1,:)' - src_pts(1,:)).^2 + (src_pts(2,:)' - src_pts(2,:)).^2);
dists = dists + diag(Inf*ones(n_src,1));
min_dists = min(dists);
disp("Min distance between points: " + num2str(min(min_dists)));


%% Do the rest.


zk = 10;
kern_0 = @(s,t) helm2d_kernel(zk, s,t);
src_info = struct;
% Source and target points are random in [-0.5, 0.5] x [-0.5, 0.5]
src_info.r = src_pts;
targ_info = struct;
targ_info.r = target_pts;

% Source weights are random uniform in [0, 1]
mu = rand(n_src, 1);
K_exact = kern_0(src_info, targ_info);

target_vals = K_exact * mu;

n_nbr = 100;

% Loop through different tolerances and record errors and timings
tol_vals = [1e-02 1e-03 1e-04 1e-05 1e-06 1e-07 1e-08 1e-09];
n_tol_vals = size(tol_vals, 2);
error_vals = zeros(n_tol_vals, 1);
dx_vals = zeros(n_tol_vals, 1);
nbinpts_vals = zeros(n_tol_vals, 1);
nproxy_vals = zeros(n_tol_vals, 1);

for i = 1:n_tol_vals
    tol = tol_vals(i);
    [grid_info, proxy_info] = get_grid(kern_0, src_info, targ_info, tol, n_nbr);
    [A_spread_s, sort_info_s ]= get_spread(kern_0, [], src_info, ...
        grid_info, proxy_info);
    [A_spread_t, sort_info_t ]= get_spread(kern_0, [], targ_info, ...
        grid_info, proxy_info);

    A_addsub = get_addsub(kern_0, [], ...
        grid_info, proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);

    k0hat = get_kernhat(kern_0,grid_info);
    evals_approx = pcfft_apply(mu,A_spread_s,A_spread_t,A_addsub,k0hat);

    % Compute relative L infinity error
    diffs = abs(evals_approx - target_vals);
    rel_linf_error = max(diffs) / max(abs(target_vals));

    % Save error and grid info.
    error_vals(i) = rel_linf_error;
    dx_vals(i) = grid_info.dx;
    nbinpts_vals(i) = grid_info.nbinpts;
    nproxy_vals(i) = proxy_info.nproxy;

    disp("tol: " + num2str(tol) + ", rel linf error: " + num2str(rel_linf_error));

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
plot(tol_vals(:), nproxy_vals(:), '.-');
hold on;
plot(tol_vals(:), nbinpts_vals(:), '.-');
% plot(tol_vals(:), nspread_vals(:), '.-');
legend("nproxy", "nbinpts");
grid on;
xscale('log');