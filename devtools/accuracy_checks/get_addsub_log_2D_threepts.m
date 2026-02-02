addpath(genpath("../../pcfft"));
close all;
clear;


rad = 2.0;


% Set up two source points
rng(4);
n_src = 2;
source_pts = [-0.0531 0.00223
                -0.00497 0.4143];

% target points are close but not exactly = source points
target_pts = [-0.0491 0.00234 0.1
                -0.00678 0.4102 0.1];
disp(size(target_pts));


src_weights = zeros(n_src,1);
src_weights(1) = 1.0;
src_weights = src_weights(:);

% Define the kernel
k = @(s,t) log_kernel(s,t);

K_src_to_target = log_kernel(struct('r',source_pts), struct('r',target_pts));

target_vals = K_src_to_target * src_weights;
n_nbr = 3; % 10000 points / 500 is approximately 20 boxes

src_info = struct;
src_info.r = source_pts;
targ_info = struct;
targ_info.r = target_pts;

%% Part 1: inspect things for 2 pts




tol = 1e-08;

[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);

[A_spread_s, K_src_to_reg, sort_info_s] = get_spread(k, k, src_info, ...
    grid_info, proxy_info);

[A_spread_t, K_targ_to_reg, sort_info_t] = get_spread(k, k, targ_info, ...
    grid_info, proxy_info);


[A_add, A_sub] = get_addsub(k, k, k, k, src_info, targ_info, grid_info, ...
    proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);

% A_addsub = A_add - A_sub;

K_grid2grid = log_kernel(grid_info, grid_info);

term1 = A_add * src_weights;
disp("main: term1: ")
disp(term1)

term2 = A_sub * src_weights;
disp("main: term2: ")
disp(term2)

AKA = A_spread_t.' * K_grid2grid * A_spread_s;

term3 = AKA * src_weights;
disp("main: term3: ")
disp(term3)

evals_approx = term1 - term2 + term3;

disp("main: evals_approx: ")
disp(evals_approx)
disp("main: target_vals: ")
disp(target_vals)

errors_at_target = max(abs(evals_approx - target_vals)) / max(abs(target_vals));
disp("errors_at_target: " + num2str(errors_at_target));


% Print out A_add, A_sub, and AKA for debugging
disp("main: A_add: ")
disp(full(A_add))
disp("main: A_sub: ")
disp(full(A_sub))
disp("main: AKA: ")
disp(full(AKA))

assert(errors_at_target < tol);


%% Part 2: Re-do the above but with a large grid where some of the points are 
% near and some are far.

tol = 1e-05;
src_info_dummy = struct;
src_info_dummy.r = (rand(2, 100) - 0.5);

targ_info_dummy = struct;
targ_info_dummy.r = (rand(2, 100) - 0.5);
[grid_info, proxy_info] = get_grid(k, src_info_dummy, targ_info_dummy, tol, n_nbr);


[A_spread_s, K_src_to_reg, sort_info_s] = get_spread(k, k, src_info, ...
    grid_info, proxy_info);

[A_spread_t, K_targ_to_reg, sort_info_t] = get_spread(k, k, targ_info, ...
    grid_info, proxy_info);

[A_add, A_sub] = get_addsub(k, k, k, k, src_info, targ_info, grid_info, ...
    proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);


K_grid2grid = log_kernel(grid_info, grid_info);

term1 = A_add * src_weights;
disp("main: term1: ")
disp(term1)

term2 = A_sub * src_weights;
disp("main: term2: ")
disp(term2)

AKA = A_spread_t.' * K_grid2grid * A_spread_s;

term3 = AKA * src_weights;
disp("main: term3: ")
disp(term3)

evals_approx = term1 - term2 + term3;

disp("main: evals_approx: ")
disp(evals_approx)
disp("main: target_vals: ")
disp(target_vals)





errors_at_target = max(abs(evals_approx - target_vals)) / max(abs(target_vals));
disp("errors_at_target: " + num2str(errors_at_target));



% % Print out A_add, A_sub, and AKA for debugging
disp("main: A_add: ")
disp(full(A_add))
disp("main: A_sub: ")
disp(full(A_sub))
disp("main: AKA: ")
disp(full(AKA))

% Print out the matrices for debugging
assert(errors_at_target < tol);
