addpath(genpath("../../pcfft"));
close all;
clear;

% Set up many random source and target points
rng(4);
n_src = 200;
n_targ = 300;
n_nbr = 100;
kern_0 = @(s,t) log_kernel(s,t);
src_info = struct;
% Source and target points are random in [-0.5, 0.5] x [-0.5, 0.5]
src_info.r = (rand(2, n_src) - 0.5);

targ_info = struct;
targ_info.r = (rand(2, n_targ) - 0.5);

src_weights = rand(n_src, 1);
K_exact = kern_0(src_info, targ_info);
target_vals = K_exact * src_weights;



tol_vals = [1e-02 1e-03 1e-04 1e-05 1e-06 1e-07 1e-08];
n_tol_vals = size(tol_vals, 2);
error_vals = zeros(n_tol_vals, 1);
% dx_vals = zeros(n_tol_vals, 1);
% nspread_vals = zeros(n_tol_vals, 1);
nproxy_vals = zeros(n_tol_vals, 1);
nbinpts_vals = zeros(n_tol_vals, 1);



for i = 1:n_tol_vals
    tol = tol_vals(i);

    [grid_info, proxy_info] = get_grid(kern_0, src_info, targ_info, tol, n_nbr);


    [A_spread_s, K_src_to_reg_s, sort_info_s ]= get_spread(kern_0, kern_0, src_info, ...
    grid_info, proxy_info);
    [A_spread_t, K_src_to_reg_t, sort_info_t ]= get_spread(kern_0, kern_0, targ_info, ...
    grid_info, proxy_info);


    [A_addsub] = get_addsub(kern_0, kern_0, kern_0, kern_0, src_info, targ_info, ...
    grid_info, proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);


    K_grid2grid = log_kernel(grid_info, grid_info);

    term1 = A_addsub * src_weights;
    term3 = A_spread_t.' * K_grid2grid * A_spread_s * src_weights;
    evals_approx = term1  + term3;
    disp("tol: " + num2str(tol) + ", evals_approx: ");
    disp(evals_approx);
    disp("target_vals: ");
    disp(target_vals);

    errors_at_target = max(abs(evals_approx(:) - target_vals(:))) / max(abs(target_vals(:)));
    disp("tol: " + num2str(tol) + ", errors_at_target: " + num2str(errors_at_target));

    % Save the error and dx vals
    error_vals(i) = errors_at_target;
    % dx_vals(i) = grid_info.dx;
    % nspread_vals(i) = grid_info.nspread;
    nproxy_vals(i) = proxy_info.nproxy;
    nbinpts_vals(i) = grid_info.nbinpts;
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
xlabel("Tolerance");
grid on;
xscale('log');
