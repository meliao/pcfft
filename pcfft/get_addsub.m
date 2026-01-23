function A_add_sub = get_addsub(kern_0, kern_s, kern_t, kern_st, src_info, targ_info, grid_info, proxy_info, sort_info_s, sort_info_t, K_src_to_reg_s, K_src_to_reg_t)


    N_bins = grid_info.nbin(1) * grid_info.nbin(2);

    id_start = src_info.id_start;
    r_sorted = src_info.r_sorted;
    % Loop through all of the bins
    for i = 1:size(id_start, 2) -1
        bin_idx = i-1; % Because bins are 0-indexed
        idx_start = id_start(i);
        idx_end = id_start(i+1) - 1;

        src_pts_in_i = r_sorted(:, idx_start:idx_end);

        % Get the regular grid points and centers of bin i
        if dim == 2
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_2d(bin_idx, grid_info);
        else
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_3d();
        end

        % Find the intersecting bins
        intersecting_bin_idxes = intersecting_bins_2d(bin_idx, grid_info, proxy_info);

        % Compute the exact near-field interactions


    end

end