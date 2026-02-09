addpath(genpath("../../pcfft"));
% close all;
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

%% Part 1: inspect things for 2 pts.
% In this test, we want to evaluate the accuracy of the addsub matrices by 
% evaluating equation (5) in the notes. On one side of this equation, we have 
% the exact interactions between sources and targets, which we can compute
% directly. On the other side, we have B mu + AKA mu.


% In this test, there are 2 source points and 3 target points. 
% Here are the dense interactions:
% s(1) -> t(1) : both in box 0
% s(2) -> t(2) : both in box 5
% s(1) -> t(3) : s(1) in box 0, t(3) in box 7. These are near.

% So here are the interactions which are NOT dense:
% s(1) -> t(2) % 0 and 5 are not near
% s(2) -> t(1) % 5 and 0 are not near
% s(2) -> t(3) % 5 and 7 are not near



tol = 1e-08;

[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);

[A_spread_s, K_src_to_reg, sort_info_s] = get_spread(k, k, src_info, ...
    grid_info, proxy_info);

% [A_spread_t, K_targ_to_reg, sort_info_t] = get_spread(k, k, targ_info, ...
%     grid_info, proxy_info);
% disp("test_neighbor_template_2d_1: Checking bin idx " + int2str(0));
% [pts_i, ~, row_idxes_i] = grid_pts_for_box_2d(0, grid_info);
% disp("row_idxes_i: ")
% disp(row_idxes_i);

% Plot the grid points with blue circles and the box grid points with red x's.
% figure;
% plot(grid_info.r(1, :), grid_info.r(2, :), 'bo');
% hold on;
% % Plot text labels <row_idxes_i> at the box grid points.
% text(pts_i(1, :), pts_i(2, :), string(row_idxes_i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

A_spread_s = A_spread_s(:, sort_info_s.ptid_srt);
% A_spread_t = A_spread_t(:, sort_info_t.ptid_srt);

n_dummy = grid_info.nbinpts^2;
N_src = size(src_info.r, 2);
n_gridpts = grid_info.ngrid(1) * grid_info.ngrid(2);
A_spread_s = [A_spread_s; sparse(n_dummy, N_src)];
dummy_idxes = n_gridpts + 1: n_gridpts + n_dummy;

% For each bin, get the neighbor bin idxes and assert that 
% indexing the rows of A_spread_s with the nbr_grididxes gets all of the relevant entries for a given bin_idx.

for bin_idx = 0 : grid_info.nbin(1)*grid_info.nbin(2)-1
    [nbr_binids, nbr_gridpts, nbr_grididxes, ~] = neighbor_template_2d(grid_info, proxy_info, bin_idx);
    
    % Get the source points in the neighbor bins
    source_idx = [];
    for j = 1:length(nbr_binids)
        source_bin_idx = nbr_binids(j);
        % disp("test_neighbor_template_2d: For bin_idx " + int2str(bin_idx) + ", neighbor bin idx " + int2str(source_bin_idx));
        

        % If the source bin idx is invalid, skip it.
        if source_bin_idx == -1
            continue;
        end
        % Source points in bin j
        idx_sj_start = sort_info_s.id_start(source_bin_idx + 1);
        idx_sj_end = sort_info_s.id_start(source_bin_idx + 2) - 1;

        source_idx = [source_idx, idx_sj_start:idx_sj_end];
    end
    

    a_1 = A_spread_s(:, source_idx);
    a_2 = A_spread_s(nbr_grididxes, source_idx);
    disp("test_neighbor_template_2d: For bin_idx " + int2str(bin_idx) + ", norm of A_spread_s(:, source_idx): " + num2str(norm(a_1, "fro")) + ", norm of A_spread_s(nbr_grididxes, source_idx): " + num2str(norm(a_2, "fro")));

    % Assert norms are close.
    assert(norm(a_1, "fro") - norm(a_2, "fro") < 1e-12);
end