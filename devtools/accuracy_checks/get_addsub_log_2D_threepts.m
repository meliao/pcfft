addpath(genpath("../../pcfft"));
close all;
clear;


rad = 2.0;


% Set up two source points
rng(4);
n_src = 2;
source_pts = [-0.0531 0.00123
                0.00497 -1.6143];

% target points are close but not exactly = source points
target_pts = [-0.0491 0.00234
                0.00678 -1.6102];


src_weights = zeros(n_src,1);
src_weights(2) = 1.0;
src_weights = src_weights(:);

% Define the kernel
k = @(s,t) log_kernel(s,t);

K_src_to_target = log_kernel(struct('r',source_pts), struct('r',target_pts));

target_vals = K_src_to_target * src_weights;
n_nbr = 10; % 10000 points / 500 is approximately 20 boxes


%% Part 1: inspect things for 2 pts

src_info = struct;
src_info.r = source_pts;

% Make the "target points" the same as the source points.
targ_info = struct;
targ_info.r = target_pts;





tol = 1e-04;

[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);

[A_spread, K_src_to_reg, sort_info_s] = get_spread(k, k, src_info, ...
    grid_info, proxy_info);

sort_info_t = SortInfo(targ_info.r, grid_info.dx, grid_info.Lbd, ...
    grid_info.nbin, grid_info.nbinpts);

[A_add, A_sub] = get_addsub(k, k, k, k, src_info, targ_info, grid_info, ...
    proxy_info, sort_info_s, sort_info_t, K_src_to_reg);

spread_weights = A_spread * src_weights;

K_reg_to_target = log_kernel(grid_info, targ_info);

term1 = A_add * src_weights;
disp("main: term1: ")
disp(term1)
term2 = A_sub * spread_weights;
disp("main: term2: ")
disp(term2)



term3 = K_reg_to_target * (A_spread * src_weights);
disp("main: term3: ")
disp(term3)

evals_approx = term1 - term2 + term3;

disp("main: evals_approx: ")
disp(evals_approx)
disp("main: target_vals: ")
disp(target_vals)





errors_at_target = max(abs(evals_approx - target_vals));
disp("errors_at_target: " + num2str(errors_at_target));


disp("main: A_add: ")
disp(full(A_add))

sub_spread = A_sub * A_spread;
disp("main: sub_spread: ")
disp(full(sub_spread))

KA = K_reg_to_target * A_spread;
disp("main: KA: ")
disp(KA)

disp("main: sub_spread - KA: ")
disp(full(sub_spread - KA))

% disp("grid_info.ngrid: ")
% disp(grid_info.ngrid)
% disp("grid_info.dx: ")
% disp(grid_info.dx)
% disp("grid_info.Lbd:")
% disp(grid_info.Lbd)

% disp("grid_info.r shape:")
% disp(size(grid_info.r))
% disp(grid_info.r(:, 1:100))

% Plot the source points and the grid points
% figure(1);
% scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
% hold on;
% scatter(src_info.r(1,:), src_info.r(2,:), 'ro');


% figure(2);
% plot(target_vals, 'bo-');
% hold on;
% plot(approx_targ_vals, 'rx-');
% legend('Exact Target Values', 'Approx Target Values');
% figure(3);
% plot(target_vals - approx_targ_vals, 'k.-');
% title('Errors at Target Points');

%% Part 2: Re-do the above but loop over tol values and plot tol vs error

% tol_vals = [1e-02 1e-03 1e-04 1e-05 1e-06 1e-07 1e-08];
% n_tol_vals = size(tol_vals, 2);
% error_vals = zeros(n_tol_vals, 1);
% dx_vals = zeros(n_tol_vals, 1);
% nspread_vals = zeros(n_tol_vals, 1);

% for i = 1:n_tol_vals
%     tol = tol_vals(i);
%     src_info = struct;
%     src_info.r = source_pts;

%     % Make the "target points" the same as the source points.
%     targ_info = struct;
%     targ_info.r = source_pts;
%     [grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);


%     % Loop through the bins and make sure that the bin center is > proxy_info.radius
%     % away from the target_pts.
%     n_bins_total = grid_info.nbin(1) * grid_info.nbin(2) ;
%     for j = 0:n_bins_total -1
%         [pts, center, row_idxes] = grid_pts_for_bin_2d(j, grid_info);
%         xdists = center(1) - target_pts(1);
%         ydists = center(2) - target_pts(2);
%         dists = sqrt(xdists.^2 + ydists.^2);
%         assert(all(dists > proxy_info.radius));
%     end
    

%     A_spread = get_spread(k, k, src_info, grid_info, proxy_info);
%     reg_weights = A_spread * src_weights;
%     reg_weights = full(reg_weights);

%     K_reg_to_target = log_kernel(grid_info, struct('r',target_pts));
%     % disp("K_reg_to_target shape: ")
%     % disp(size(K_reg_to_target))
%     % disp("reg_weights shape:")
%     % disp(size(reg_weights))

%     approx_targ_vals = K_reg_to_target * reg_weights;

%     errors_at_target = max(abs(approx_targ_vals - target_vals));

%     % Save the error and dx vals
%     error_vals(i) = errors_at_target;
%     dx_vals(i) = grid_info.dx;
%     nspread_vals(i) = grid_info.nspread;
% end


% figure(1);
% plot(tol_vals, error_vals, 'o-');
% hold on;
% plot(tol_vals, tol_vals, 'k--');
% xscale('log');
% yscale('log');
% xlabel("Error tolerance");
% ylabel("Observed error");

% figure(2);
% plot(tol_vals, dx_vals, "o-");
% hold on;
% xscale('log');
% xlabel("Error tolerance");
% ylabel("dx");

% figure(3);
% plot(tol_vals, nspread_vals, "o-");
% hold on;
% xscale('log');
% xlabel("Error tolerance");
% ylabel("nspread");