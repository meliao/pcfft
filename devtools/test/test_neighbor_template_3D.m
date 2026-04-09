% Makes sure get_addsub returns without error on a 2D input.
addpath(genpath('../../pcfft'));

close all;

tol = 1e-6;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 1037;
rng(0);
src_info_2d.r = (rand(3, n_src) - 0.5);
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
ntarg = 173;
targ_info_2d.r = rand(3, ntarg) - 0.5;


[grid_info, proxy_info] = get_grid(@one_over_r_kernel, ...
    src_info_2d, targ_info_2d, tol);

N_bin = grid_info.nbin(1) * grid_info.nbin(2);

% disp("test_intersecting_bins_2d: N_bin = " + int2str(N_bin));

[~, tmpl_pts, tmpl_idxes] = abstract_neighbor_spreading_3D(grid_info, proxy_info);

disp("test: tmpl_pts size: " + int2str(size(tmpl_pts)));
disp("test: tmpl_idxes size: " + int2str(size(tmpl_idxes)));
[nbr_binids, nbr_gridpts, nbr_grididxes] = neighbor_template_3d(grid_info, proxy_info, 8, tmpl_pts, tmpl_idxes);


% nbr_grididxes should either be in the range [1, ngrid(1)*ngrid(2)] or be equal to dummy_idx = ngrid(1)*ngrid(2) + 1
dummy_idx = grid_info.ngrid(1) * grid_info.ngrid(2) * grid_info.ngrid(3) + 1;
assert(all(nbr_grididxes >= 1 & nbr_grididxes <= dummy_idx));

% valid nbr_grididxes should correctly index nbr_gridpts
keep_bool = nbr_grididxes ~= dummy_idx;
valid_grididxes = nbr_grididxes(keep_bool);
valid_gridpts = nbr_gridpts(:, keep_bool);


% Scatter plot the valid nbr_gridpts
scatter3(grid_info.r(1,:), grid_info.r(2,:), grid_info.r(3,:), 10, 'b');
hold on;
scatter3(valid_gridpts(1,:), valid_gridpts(2,:), valid_gridpts(3,:), 20, 'r', 'filled');

% Plot the text of nbr_binids at the center of each bin
for i = 1:size(nbr_binids, 2)
    bin_id = nbr_binids(i);
    ctr = bin_center(bin_id, grid_info);
    text(ctr(1), ctr(2), ctr(3), string(bin_id), 'Color', 'k');
end

title("Grid points with their linear indices");

N_bins = grid_info.nbin(1) * grid_info.nbin(2) * grid_info.nbin(3);
for bin_idx = 0:N_bins-1
    [nbr_binids, nbr_gridpts, nbr_grididxes] = neighbor_template_3d(grid_info, proxy_info, bin_idx, tmpl_pts, tmpl_idxes);
    assert(all(nbr_grididxes >= 1 & nbr_grididxes <= dummy_idx));

    keep_bool = nbr_grididxes ~= dummy_idx;
    valid_grididxes = nbr_grididxes(keep_bool);

    % Assert valid_grididxes are unique
    assert(length(unique(valid_grididxes)) == length(valid_grididxes));

    % Assert valid_grididxes correctly index valid_gridpts
    valid_gridpts = nbr_gridpts(:, keep_bool);
    for i = 1:size(valid_grididxes, 2)
        idx = valid_grididxes(i);
        pt = valid_gridpts(:, i);
        grid_pt = grid_info.r(:, idx);
        dist = norm(pt - grid_pt);
        assert(dist < 1e-12);
    end
end

close all;