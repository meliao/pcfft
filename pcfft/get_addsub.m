function [A_add, A_sub] = get_addsub(kern_0, kern_s, kern_t, kern_st, src_info, ...
    targ_info, grid_info, proxy_info, sort_info_s, sort_info_t, ...
    A_spread_s, A_spread_t)

    N_src = size(src_info.r, 2);
    N_targ = size(targ_info.r, 2);

    max_bin_idx = grid_info.nbin(1) * grid_info.nbin(2) - 1;
    n_dummy = grid_info.nbinpts^2;
    n_gridpts = grid_info.ngrid(1) * grid_info.ngrid(2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build a spreading template matrix for adjacent source points.
    % Then build a list of regular gridpoints that are in the intersecting bins
    [reg_neighbor_template_pts, nbr_bin_idx] = neighbor_template_2d(grid_info, proxy_info);
    nbr_info = struct('r', reg_neighbor_template_pts);
    [pts0, ctr_0, ~] = grid_pts_for_bin_2d(nbr_bin_idx, grid_info);
    % pts0_centered = pts0 - ctr_0;
    bin_info = struct('r', pts0);
    K_nbr2bin = kern_0(nbr_info, bin_info);

    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.
    A_add = sparse(N_targ, N_src);
    A_sub = sparse(N_targ, N_src);

    % Sort the cols of A_spread_s to match the sorted source points
    A_spread_s = A_spread_s(:, sort_info_s.ptid_srt);

    % Add n_dummy rows of zeros to A_spread_s to handle empty bins
    A_spread_s = [A_spread_s; sparse(n_dummy, N_src)];
    dummy_idxes = n_gridpts + 1: n_gridpts + n_dummy;
    % Sort the cols of A_spread_t to match the sorted target points
    A_spread_t = A_spread_t(:, sort_info_t.ptid_srt);

    % Loop through all of the bins. 
    for i = 1:size(sort_info_s.id_start, 2) -1
        bin_idx = i - 1; % Because bins are 0-indexed
        disp("get_addsub: Processing bin " + int2str(bin_idx));

        % Need the center of bin i to center the source points, and need the 
        % indexes of the regular grid points for spreading bin i, so we can
        % correctly index A_spread_t.
        [~, ctr_i, reg_idxs_i] = grid_pts_for_bin_2d(bin_idx, grid_info);

        % Target points in bin i
        idx_ti_start = sort_info_t.id_start(i);
        idx_ti_end = sort_info_t.id_start(i + 1) - 1;
        
        targ_pts_in_i = sort_info_t.r_srt(:, idx_ti_start:idx_ti_end);

        % Find the intersecting bins
        [id_xs, id_ys, nbr_binids] = intersecting_bins_2d(bin_idx, grid_info, proxy_info);
        % disp("get_addsub: nbr_binids: ");
        % disp(nbr_binids);

        % Loop through all of the neighbor bins and fill in the local source points. 
        % After this loop, we will update A_add and A_sub with the neigbors of bin i.
        source_loc = [];
        source_idx = [];
        reg_idxes = [];

        for j = 1:size(nbr_binids, 2)
            % This iter of the loop does interaction between target bin i and 
            % source bin j
            source_bin_idx = nbr_binids(j);

            % If the source bin idx is invalid, skip it.
            if source_bin_idx == -1
                reg_idxes = [reg_idxes, dummy_idxes];
                continue;
            end

            % Source points in bin j
            idx_sj_start = sort_info_s.id_start(source_bin_idx + 1);
            idx_sj_end = sort_info_s.id_start(source_bin_idx + 2) - 1;

            loc_src_in_j = sort_info_s.r_srt(:, idx_sj_start:idx_sj_end);

            source_loc = [source_loc, loc_src_in_j];
            source_idx = [source_idx, idx_sj_start:idx_sj_end];

            % Get the regular grid point idxes for spreading bin j
            [~, ~, reg_idxs_j] = grid_pts_for_bin_2d(source_bin_idx, grid_info);
            reg_idxes = [reg_idxes, reg_idxs_j];

        end

        % It may be the case that there are no source points in the bins 
        % neighboring target bin i. 
        if isempty(source_idx)
            continue;
        end

        % disp("get_addsub: bin " + int2str(bin_idx) + " source_loc size: ");
        % disp(size(source_loc));

        % Update A_add with exact near-field interactions.
        K_src2targ = kern_0(struct('r', source_loc), ...
                            struct('r', targ_pts_in_i ));
        A_add(idx_ti_start:idx_ti_end, source_idx) = ...
            A_add(idx_ti_start:idx_ti_end, source_idx) + K_src2targ;

        % Update A_sub with approximated near-field interactions.
        A_spread_t_i = A_spread_t(reg_idxs_i, idx_ti_start:idx_ti_end);
        A_spread_s_j = A_spread_s(reg_idxes, source_idx);

        % Print shape info for A_spread_t_i
        disp("get_addsub: A_spread_t_i size: ");
        disp(size(A_spread_t_i));
        disp("get_addsub: K_nbr2bin size: ");
        disp(size(K_nbr2bin));
        disp("get_addsub: A_spread_s_j size: ");
        disp(size(A_spread_s_j));

        AKA_chunk = A_spread_t_i.' * K_nbr2bin * A_spread_s_j;

        disp("get_addsub: AKA_chunk: ")
        disp(AKA_chunk);

        A_sub(idx_ti_start:idx_ti_end, source_idx) = ...
            A_sub(idx_ti_start:idx_ti_end, source_idx) + AKA_chunk;

    end

    % A_addsub = A_add - A_sub;

    % Reorder the rows to match the original target point ordering
    A_add(sort_info_t.ptid_srt, :) = A_add;
    A_sub(sort_info_t.ptid_srt, :) = A_sub;

    % % Reorder the columns to match the original source point ordering
    A_add(:, sort_info_s.ptid_srt) = A_add(:, sort_info_s.ptid_srt);
    A_sub(:, sort_info_s.ptid_srt) = A_sub(:, sort_info_s.ptid_srt);

end
