function [A_addsub] = get_addsub(kern_0, kern_s, kern_t, kern_st, src_info, ...
    targ_info, grid_info, proxy_info, sort_info_s, sort_info_t, ...
    A_spread_s, A_spread_t)

    N_src = size(src_info.r, 2);
    N_targ = size(targ_info.r, 2);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build a spreading template matrix for adjacent source points.
    % Then build a list of regular gridpoints that are in the intersecting bins
    [~, reg_neighbor_template_pts, ~, nbr_bin_idx] = neighbor_template_2d(grid_info, proxy_info);
    nbr_info = struct('r', reg_neighbor_template_pts);
    [pts0, ctr_0, ~] = grid_pts_for_box_2d(nbr_bin_idx, grid_info);
    % pts0_centered = pts0 - ctr_0;
    bin_info = struct('r', pts0);
    K_nbr2bin = kern_0(nbr_info, bin_info);

    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.
    % A_addsub = sparse(N_targ, N_src);

    % Sort the cols of A_spread_s and A_spread_t to match the sorted source points
    A_spread_s = A_spread_s(:, sort_info_s.ptid_srt);
    A_spread_t = A_spread_t(:, sort_info_t.ptid_srt);


    % Add 1 row of zeros to A_spread_s to handle empty bins
    A_spread_s = [A_spread_s; sparse(1, N_src)];
    % dummy_idx = n_gridpts + 1;

    % TODO: correct formula for number of corrections
    ncor = grid_info.n_nbr*ceil(mean([N_targ,N_src]));

    % These are the arrays we will use to build the sparse A_addsub
    % in COO format.
    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.
    iid = zeros(1,ncor);
    jid = zeros(1,ncor);
    vals = zeros(1,ncor);
    id_start = 0;


    % Loop through all of the bins. 
    for i = 1:size(sort_info_s.id_start, 2) -1
        bin_idx = i - 1; % Because bins are 0-indexed
        % disp("get_addsub: Processing bin " + int2str(bin_idx));

        % Need the center of bin i to center the source points, and need the 
        % indexes of the regular grid points for spreading bin i, so we can
        % correctly index A_spread_t.
        [~, ~, reg_idxs_i] = grid_pts_for_box_2d(bin_idx, grid_info);

        % Target points in bin i
        idx_ti_start = sort_info_t.id_start(i);
        idx_ti_end = sort_info_t.id_start(i + 1) - 1;
        
        targ_pts_in_i = sort_info_t.r_srt(:, idx_ti_start:idx_ti_end);

        % Build the spreading template
        [nbr_binids, ~, nbr_grididxes, ~] = ...
            neighbor_template_2d(grid_info, proxy_info, bin_idx);

        % Loop through all of the neighbor bins and fill in the local source points. 
        % After this loop, we will update A_add and A_sub with the neigbors of bin i.
        % TODO: preallocate this array. We can infer the length from sort_info_s.
        source_loc = [];
        source_idx = [];

        for j = 1:length(nbr_binids)
            % This iter of the loop does interaction between target bin i and 
            % source bin j
            source_bin_idx = nbr_binids(j);

            % If the source bin idx is invalid, skip it.
            if source_bin_idx == -1
                continue;
            end

            % Source points in bin j
            idx_sj_start = sort_info_s.id_start(source_bin_idx + 1);
            idx_sj_end = sort_info_s.id_start(source_bin_idx + 2) - 1;

            loc_src_in_j = sort_info_s.r_srt(:, idx_sj_start:idx_sj_end);

            source_loc = [source_loc, loc_src_in_j];
            source_idx = [source_idx, idx_sj_start:idx_sj_end];


        end

        % It may be the case that there are no source points in the bins 
        % neighboring target bin i. 
        if isempty(source_idx)
            continue;
        end

        % Update A_addsub with exact near-field interactions. This is the "add"
        % part.
        K_src2targ = kern_0(struct('r', source_loc), ...
                            struct('r', targ_pts_in_i ));

        % Update A_sub with approximated near-field interactions. This is the 
        % "sub" part.
        A_spread_t_i = A_spread_t(reg_idxs_i, idx_ti_start:idx_ti_end);
        A_spread_s_j = A_spread_s(nbr_grididxes, source_idx);
        AKA_chunk = A_spread_t_i.' * K_nbr2bin * A_spread_s_j;

        Aloc =  K_src2targ - AKA_chunk;

        % A_addsub(idx_ti_start:idx_ti_end, source_idx) = Aloc;

        % Update COO arrays.
        is = (idx_ti_start:idx_ti_end);
        js = source_idx;
        is = repmat(is(:), 1, size(js,2));
        js = repmat(js(:).', size(is,1), 1);
        n_sparse = numel(Aloc);

        iid(id_start + (1:n_sparse)) = is(:).';
        jid(id_start + (1:n_sparse)) = js(:).';
        vals(id_start + (1:n_sparse)) = Aloc(:).';
        id_start = id_start + n_sparse;

    end

    A_addsub = sparse(iid(1:id_start), jid(1:id_start), vals(1:id_start), N_targ, N_src);

    % Reorder the rows to match the original target point ordering
    A_addsub(sort_info_t.ptid_srt, :) = A_addsub;

    % % Reorder the columns to match the original source point ordering
    A_addsub(:, sort_info_s.ptid_srt) = A_addsub; 
end
