addpath(genpath("../../pcfft"));
close all;
clear;

% Set up random source and target points
rng(1);
n_src = 500;
n_targ = 517;
dim = 3;

kern_0 = @(s,t) one_over_r_kernel(s,t);
src_info = struct;
% Source and target points are random in [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]
src_info.r = (rand(dim, n_src) - 0.5);
targ_info = struct;
targ_info.r = (rand(dim, n_targ) - 0.5);

% Source weights are random uniform in [0, 1]
mu = rand(n_src, 1);
K_exact = kern_0(src_info, targ_info);
target_vals = K_exact * mu;

tol = 1e-6;
n_nbr = 100;
[grid_info, proxy_info] = get_grid(kern_0, src_info, targ_info, tol, n_nbr);

[A_spread_s, sort_info_s ]= get_spread(kern_0, [], src_info, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t ]= get_spread(kern_0, [], targ_info, ...
    grid_info, proxy_info);

A_addsub = get_addsub(kern_0, [], grid_info, proxy_info, sort_info_s, ...
    sort_info_t, A_spread_s, A_spread_t);

k0hat = get_kernhat(kern_0,grid_info);
evals_approx = pcfft_apply(mu,A_spread_s,A_spread_t,A_addsub,k0hat);

% Compute relative L infinity error
diffs = abs(evals_approx - target_vals);
rel_linf_error = max(diffs) / max(abs(target_vals));

assert(rel_linf_error < tol);