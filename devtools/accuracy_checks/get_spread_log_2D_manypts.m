addpath(genpath("../../pcfft"));

target_rad = 2.5;
ntarg = 200;
target_pts = get_ring_points(ntarg, target_rad);



% Set up two source points
rng(4);
n_src = 200;
rng(0);
source_pts = (rand(2, n_src) - 0.5) * 2;


src_weights = randn(n_src,1);
src_weights = src_weights(:);

% Define the kernel
k = @(s,t) log_kernel(s,t);

K_src_to_target = log_kernel(source_pts, target_pts);

target_vals = K_src_to_target * src_weights;
n_nbr = 500; % 10000 points / 500 is approximately 20 boxes


%% Part 1: inspect things for 2 pts

% src_info = struct;
% src_info.r = source_pts;

% % Make the "target points" the same as the source points.
% targ_info = struct;
% targ_info.r = source_pts;





% tol = 1e-07;

% [grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);


% % Beef up proxy_info
% % nproxy = proxy_info.n_points_total * 10;
% % proxy_pts = get_ring_points(nproxy, proxy_info.radius);
% % proxy_info.r = proxy_pts;
% % proxy_info.n_points_total = nproxy;

% nbin = grid_info.nspread - 1;

% % disp("grid_info.ngrid: ")
% % disp(grid_info.ngrid)
% % disp("grid_info.dx: ")
% % disp(grid_info.dx)
% % disp("grid_info.Lbd:")
% % disp(grid_info.Lbd)

% % disp("grid_info.r shape:")
% % disp(size(grid_info.r))
% % disp(grid_info.r(:, 1:100))

% % Plot the source points and the grid points
% figure(1);
% scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
% hold on;
% scatter(src_info.r(1,:), src_info.r(2,:), 'ro');

% A_spread = get_spread(k, k, src_info, grid_info, proxy_info, nbin);

% % disp("A_spread shape:")
% % disp(size(A_spread))

% reg_weights = A_spread * src_weights;
% reg_weights = full(reg_weights);

% K_reg_to_target = log_kernel(grid_info.r, target_pts);
% % disp("K_reg_to_target shape: ")
% % disp(size(K_reg_to_target))
% % disp("reg_weights shape:")
% % disp(size(reg_weights))

% approx_targ_vals = K_reg_to_target * reg_weights;

% errors_at_target = max(abs(approx_targ_vals - target_vals));
% % disp("target_vals: " + num2str(target_vals'));
% % disp("approx_targ_vals: " + num2str(approx_targ_vals'));
% disp("errors_at_target: " + num2str(errors_at_target));

% figure(2);
% plot(target_vals, 'bo-');
% hold on;
% plot(approx_targ_vals, 'rx-');
% legend('Exact Target Values', 'Approx Target Values');
% figure(3);
% plot(target_vals - approx_targ_vals, 'k.-');
% title('Errors at Target Points');

%% Part 2: Re-do the above but loop over tol values and plot tol vs error

tol_vals = [1e-02 1e-03 1e-04 1e-05 1e-06 1e-07 1e-08];
n_tol_vals = size(tol_vals, 2);
error_vals = zeros(n_tol_vals, 1);
dx_vals = zeros(n_tol_vals, 1);
nspread_vals = zeros(n_tol_vals, 1);

for i = 1:n_tol_vals
    tol = tol_vals(i);
    src_info = struct;
    src_info.r = source_pts;

    % Make the "target points" the same as the source points.
    targ_info = struct;
    targ_info.r = source_pts;
    [grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);
    nbin = grid_info.nspread -1;


    A_spread = get_spread(k, k, src_info, grid_info, proxy_info, nbin);
    reg_weights = A_spread * src_weights;
    reg_weights = full(reg_weights);

    K_reg_to_target = log_kernel(grid_info.r, target_pts);
    % disp("K_reg_to_target shape: ")
    % disp(size(K_reg_to_target))
    % disp("reg_weights shape:")
    % disp(size(reg_weights))

    approx_targ_vals = K_reg_to_target * reg_weights;

    errors_at_target = max(abs(approx_targ_vals - target_vals));

    % Save the error and dx vals
    error_vals(i) = errors_at_target;
    dx_vals(i) = grid_info.dx;
    nspread_vals(i) = grid_info.nspread;
end


figure(1);
plot(tol_vals, error_vals, 'o-');
hold on;
plot(tol_vals, tol_vals, 'k--');
xscale('log');
yscale('log');
xlabel("Error tolerance");
ylabel("Observed error");

figure(2);
plot(tol_vals, dx_vals, "o-");
hold on;
xscale('log');
xlabel("Error tolerance");
ylabel("dx");

figure(3);
plot(tol_vals, nspread_vals, "o-");
hold on;
xscale('log');
xlabel("Error tolerance");
ylabel("nspread");

figure(4);
scatter(source_pts(1,:), source_pts(2,:));
