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
n_nbr = 5;

[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);

[A_spread_s, sort_info_s] = get_spread(k, k, src_info, ...
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

[~, tmpl_pts, tmpl_idxes] = abstract_neighbor_spreading_2D(grid_info, proxy_info);
for bin_idx = 0 : grid_info.nbin(1)*grid_info.nbin(2)-1
    [nbr_binids, nbr_gridpts, nbr_grididxes] = neighbor_template_2d(grid_info, proxy_info, bin_idx, tmpl_pts, tmpl_idxes);

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
    % disp("test_neighbor_template_2d: For bin_idx " + int2str(bin_idx) + ", norm of A_spread_s(:, source_idx): " + num2str(norm(a_1, "fro")) + ", norm of A_spread_s(nbr_grididxes, source_idx): " + num2str(norm(a_2, "fro")));

    % Assert norms are close.
    assert(norm(a_1, "fro") - norm(a_2, "fro") < 1e-12);
end

%% Part 2: many pts

n_src = 10000;
n_targ = 5017;
dim = 2;

kern_0 = @(s,t) log_kernel(s,t);
src_info = struct;
% Source and target points are random in [-0.5, 0.5] x [-0.5, 0.5]
src_info.r = (rand(dim, n_src) - 0.5);
targ_info = struct;
targ_info.r = (rand(dim, n_targ) - 0.5);

tol = 1e-10;
n_nbr = 100;
[grid_info, proxy_info] = get_grid(kern_0, src_info, targ_info, tol, n_nbr);

% Loop through all of the boxes
[~, tmpl_pts, tmpl_idxes] = abstract_neighbor_spreading_2D(grid_info, proxy_info);
for bin_idx = 0 : grid_info.nbin(1)*grid_info.nbin(2)-1
    [nbr_binids, nbr_gridpts, nbr_grididxes] = neighbor_template_2d(grid_info, proxy_info, bin_idx, tmpl_pts, tmpl_idxes);

    valid_idxes = nbr_grididxes <= grid_info.ngrid(1) * grid_info.ngrid(2);
    valid_nbr_gridpts = nbr_gridpts(:, valid_idxes);
    valid_nbr_grididxes = nbr_grididxes(valid_idxes);

    % Assert that the valid grid points match the expected grid points
    for i = 1:size(valid_nbr_grididxes, 2)
        idx = valid_nbr_grididxes(i);
        pt = valid_nbr_gridpts(:, i);
        grid_pt = grid_info.r(:, idx);
        dist = norm(pt - grid_pt);
        assert(dist < 1e-12);
    end

    % Assert that # binids returned = # grid pts returned = # grid idxes returned
    assert(length(nbr_grididxes) == size(nbr_gridpts, 2));

    % Assert that the grid points are within the valid range unless marked
    % with a dummy idx.

end

%% Part 3: Plot neighbor template and proxy circle for a particular bin

n_src_3 = 10000;
n_targ_3 = 5017;
dim_3 = 2;

kern_3 = @(s,t) log_kernel(s,t);
src_info_3 = struct;
src_info_3.r = (rand(dim_3, n_src_3) - 0.5);
targ_info_3 = struct;
targ_info_3.r = (rand(dim_3, n_targ_3) - 0.5);

tol_3 = 1e-8;
n_nbr_3 = 100;
[grid_info_3, proxy_info_3] = get_grid(kern_3, src_info_3, targ_info_3, tol_3, n_nbr_3);

% Use the center bin as the bin to plot
plot_bin_idx = grid_info_3.center_bin;

% Compute the neighbor template for the chosen bin
[~, tmpl_pts_3, tmpl_idxes_3] = abstract_neighbor_spreading_2D(grid_info_3, proxy_info_3);
[~, nbr_gridpts_3, nbr_grididxes_3] = neighbor_template_2d(grid_info_3, proxy_info_3, plot_bin_idx, tmpl_pts_3, tmpl_idxes_3);

% Discard out-of-bounds dummy points
valid_mask_3 = nbr_grididxes_3 <= grid_info_3.ngrid(1) * grid_info_3.ngrid(2);
valid_nbr_gridpts_3 = nbr_gridpts_3(:, valid_mask_3);

% Center of the chosen bin and proxy radius
bin_ctr_3 = bin_center(plot_bin_idx, grid_info_3);
proxy_rad_3 = proxy_info_3.radius;

% Collect all bin centers
n_bins_3 = grid_info_3.nbin(1) * grid_info_3.nbin(2);
all_bin_ctrs_3 = zeros(2, n_bins_3);
for bi = 0 : n_bins_3 - 1
    all_bin_ctrs_3(:, bi+1) = bin_center(bi, grid_info_3);
end

% Proxy circle via get_ring_points
proxy_circle_3 = get_ring_points(300, proxy_rad_3, bin_ctr_3);


% Second proxy circle inflated by the amount used in interaction_radius
inflation_3 = sqrt(grid_info_3.dim) * grid_info_3.rpad * grid_info_3.dx;
proxy_circle_inflated_3 = get_ring_points(300, proxy_rad_3 + inflation_3, bin_ctr_3);

figure;
% Regular grid points first, in black
plot(grid_info_3.r(1,:), grid_info_3.r(2,:), 'k.', 'MarkerSize', 3);
hold on;
% Neighbor template grid points in red
plot(valid_nbr_gridpts_3(1,:), valid_nbr_gridpts_3(2,:), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
% All bin centers as stars
plot(all_bin_ctrs_3(1,:), all_bin_ctrs_3(2,:), 'g*', 'MarkerSize', 6, 'LineWidth', 1);
% Chosen bin center highlighted
plot(bin_ctr_3(1), bin_ctr_3(2), 'b*', 'MarkerSize', 12, 'LineWidth', 2);
% Proxy circle
plot([proxy_circle_3(1,:), proxy_circle_3(1,1)], [proxy_circle_3(2,:), proxy_circle_3(2,1)], 'b-', 'LineWidth', 2);
plot([proxy_circle_inflated_3(1,:), proxy_circle_inflated_3(1,1)], [proxy_circle_inflated_3(2,:), proxy_circle_inflated_3(2,1)], 'b--', 'LineWidth', 2);
hold off;
axis equal;
title(sprintf('Neighbor template for bin %d', plot_bin_idx));
legend('Regular grid pts', 'Neighbor template pts', 'All bin centers', 'Chosen bin center', 'Proxy circle', 'Location', 'best');
close all;