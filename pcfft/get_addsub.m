function [A_add, A_sub] = get_addsub(kern_0, kern_s, kern_t, kern_st, src_info, ...
    targ_info, grid_info, proxy_info, sort_info_s, sort_info_t, ...
    K_src_to_reg_s)

    N_src = size(src_info.r, 2);
    N_targ = size(targ_info.r, 2);

    N_bins = grid_info.nbin(1) * grid_info.nbin(2);

    N_reg = size(grid_info.r, 2);


    % Rows of A_add and A_sub are ordered according to sorted target points.
    % Cols of A_add and A_sub are ordered according to sorted source points.
    A_add = sparse(N_targ, N_src);
    A_sub = sparse(N_targ, N_reg);

    % Loop through all of the bins
    for i = 1:size(sort_info_s.id_start, 2) -1
        bin_idx = i - 1; % Because bins are 0-indexed

        idx_start_t = sort_info_t.id_start(i);
        idx_end_t = sort_info_t.id_start(i + 1) - 1;
        
        % Print out the idx_start_t and idx_end_t
        disp("get_addsub: Processing target bin " + int2str(bin_idx) + ...
            ", idx_start_t: " + int2str(idx_start_t) + ", idx_end_t: " + int2str(idx_end_t));

        targ_pts_in_i = sort_info_t.r_srt(:, idx_start_t:idx_end_t);

        % % Get the regular grid points and centers of bin i
        % if grid_info.dim == 2
        %     [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_2d(bin_idx, grid_info);
        % else
        %     [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_3d();
        % end

        % Find the intersecting bins
        intersecting_bin_idxes = intersecting_bins_2d(bin_idx, grid_info, proxy_info);

        % disp("get_addsub:   Intersecting bins: ");
        disp(intersecting_bin_idxes);
        % Compute the exact near-field interactions (This is for A_add)
        for j = 1:size(intersecting_bin_idxes, 2)
            % This iter of the loop does interaction between target bin i and 
            % source bin j
            source_bin_idx = intersecting_bin_idxes(j);
            % Get the source points in this bin
            idx_start_s = sort_info_s.id_start(source_bin_idx + 1);
            idx_end_s = sort_info_s.id_start(source_bin_idx + 2) - 1;

            % disp("get_addsub:   Processing source bin " + int2str(source_bin_idx) + ...
            %     ", idx_start_s: " + int2str(idx_start_s) + ", idx_end_s: " + int2str(idx_end_s));
            src_pts_in_j = sort_info_s.r_srt(:, idx_start_s:idx_end_s);

            K_src_to_targ = kern_0(struct('r', src_pts_in_j), ...
                                struct('r', targ_pts_in_i));

            A_add(idx_start_t:idx_end_t, idx_start_s:idx_end_s) = ...
                A_add(idx_start_t:idx_end_t, idx_start_s:idx_end_s) + K_src_to_targ;


            % Compute the near-field interactions between targets and 
            % regular grid points (This is for A_sub)

            if grid_info.dim == 2
                [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_2d(bin_idx, grid_info);
            else
                [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_3d();
            end

            K_reg_to_targ = kern_t(struct('r', pts_i), ...
                                struct('r', targ_pts_in_i));

            A_sub(idx_start_t:idx_end_t, row_idxes_i) = ...
                A_sub(idx_start_t:idx_end_t, row_idxes_i) + ...
                K_reg_to_targ;
        end

    end

end
