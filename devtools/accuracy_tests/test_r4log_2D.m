addpath(genpath("../../pcfft"));
close all;
clear;

% Set up random source and target points
rng(1);
n_src = 5000;
n_targ = 5017;
dim = 2;

kern_0 = @(s,t) r4log_kernel(s,t);
src_info = struct;
% Source and target points are random in [-0.5, 0.5] x [-0.5, 0.5]
src_info.r = (rand(dim, n_src) - 0.5);
targ_info = struct;
targ_info.r = (rand(dim, n_targ) - 0.5);

% Source weights are random uniform in [0, 1]
mu = rand(n_src, 1);
K_exact = kern_0(src_info, targ_info);
target_vals = K_exact * mu;

tol = 1e-10;
n_nbr = 100;
% Green's function for a sixth order PDE requires us to proxy against
% several derivatives
opts = [];
opts.proxy_der = 1;
opts.multi_shell = 0;
[grid_info, proxy_info] = get_grid(kern_0, src_info, targ_info, tol, n_nbr,opts);

[A_spread_s, sort_info_s ]= get_spread(kern_0, [], src_info, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t ]= get_spread(kern_0, [], targ_info, ...
    grid_info, proxy_info);

A_addsub = get_addsub(kern_0, [], src_info, targ_info, ...
    grid_info, proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);

k0hat = get_kernhat(kern_0,grid_info);
evals_approx = pcfft_apply(mu,A_spread_s,A_spread_t,A_addsub,k0hat);

% Compute relative L infinity error
diffs = abs(evals_approx - target_vals);
rel_linf_error = max(diffs) / max(abs(target_vals));

disp("test_r4log_2D: tol: " + num2str(tol) + ", rel linf error: " + num2str(rel_linf_error));


assert(rel_linf_error < tol);

function k_evals = r4log_kernel(src_pts,target_pts)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% Computes log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts.r(1, :) - target_pts.r(1, :).';
ry = src_pts.r(2, :) - target_pts.r(2, :).';

dist = sqrt(rx.^2 + ry.^2);

k_evals = dist.^2.*log(dist);

end