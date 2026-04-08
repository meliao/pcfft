% This tests the K_nbr2bin object constructed in get_addsub().
addpath(genpath('../../pcfft'));

k   = @(s,t) log_kernel(s,t);
tol = 1e-8;
rng(0);
src_info.r = rand(2,50) - 0.5;
[grid_info, proxy_info] = get_grid(k, src_info, src_info, tol);

[pts0, reg_neighbor_template_pts, reg_neighbor_template_idxes] = abstract_neighbor_spreading_2D(grid_info, proxy_info);
box_center = bin_center(grid_info.center_bin, grid_info);
pts0_abs = pts0; % pts0 is already absolute (before subtraction in get_addsub)

% K_nbr2bin as built in get_addsub (relative coords, center bin)
pts0_rel = pts0 - box_center;
K_nbr2bin = k(struct('r', reg_neighbor_template_pts), struct('r', pts0_rel));
K_nbr2bin_r = 0;
for d = 1:2
    K_nbr2bin_r = K_nbr2bin_r + (reg_neighbor_template_pts(d,:) - pts0_rel(d,:).').^2;
end
K_nbr2bin(K_nbr2bin_r < 1e-14) = 0;

% For each non-boundary bin, compute K directly from absolute grid coords
% and compare to K_nbr2bin
for bin_idx = 0 : grid_info.nbin(1)*grid_info.nbin(2) - 1
    [~, ~, reg_idxs_i]    = grid_pts_for_box_2d(bin_idx, grid_info);
    [~, nbr_gridpts, nbr_grididxes] = neighbor_template_2d(grid_info, proxy_info, bin_idx, reg_neighbor_template_pts, reg_neighbor_template_idxes);

    % Only consider in-bounds template points
    dummy = grid_info.ngrid(1)*grid_info.ngrid(2) + 1;
    valid = nbr_grididxes ~= dummy;

    % Absolute grid coordinates
    box_abs   = grid_info.r(:, reg_idxs_i);          % [2 x nbox]
    tmpl_abs  = grid_info.r(:, nbr_grididxes(valid)); % [2 x nvalid]

    % Direct K (absolute coords)
    K_direct = k(struct('r', tmpl_abs), struct('r', box_abs));
    r2 = 0;
    for d = 1:2
        r2 = r2 + (tmpl_abs(d,:) - box_abs(d,:).').^2;
    end
    K_direct(r2 < 1e-14) = 0;

    % K_nbr2bin restricted to valid columns
    K_approx = K_nbr2bin(:, valid);

    err = max(abs(K_direct(:) - K_approx(:)));
    assert(err < 1e-10, sprintf('bin %d: K ordering mismatch, err=%g', bin_idx, err));
end

%% test_2

bin_idx = grid_info.center_bin;
[~, tmpl_pts_abs_1, tmpl_idxes] = abstract_neighbor_spreading_2D(grid_info, proxy_info);
box_ctr = bin_center(grid_info.center_bin, grid_info);
tmpl_pts_abs = tmpl_pts_abs_1 + box_ctr; % shift from relative to absolute

[~, ~, nbr_grididxes_c] = neighbor_template_2d(grid_info, proxy_info, bin_idx, tmpl_pts_abs_1, tmpl_idxes);

% For valid points, check that nbr_grididxes index the same physical points
% as tmpl_idxes
dummy = grid_info.ngrid(1)*grid_info.ngrid(2) + 1;
valid = nbr_grididxes_c ~= dummy;
in_bounds = tmpl_idxes(1,:) >= 1 & tmpl_idxes(1,:) <= grid_info.ngrid(1) & ...
            tmpl_idxes(2,:) >= 1 & tmpl_idxes(2,:) <= grid_info.ngrid(2);

linear_from_tmpl = (tmpl_idxes(1,in_bounds)-1)*grid_info.ngrid(2) + tmpl_idxes(2,in_bounds);
assert(isequal(sort(nbr_grididxes_c(valid)), sort(linear_from_tmpl)), ...
    'center_bin: nbr_grididxes set does not match tmpl_idxes set');
assert(isequal(nbr_grididxes_c(valid), linear_from_tmpl), ...
    'center_bin: nbr_grididxes ORDER does not match tmpl_idxes order');