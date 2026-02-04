function [A_addsub] = get_addsub(kern_0, kern_st, src_info, ...
    targ_info, grid_info, proxy_info, sort_info_s, sort_info_t, ...
    A_spread_s, A_spread_t)

    der_fields_s = fieldnames(sort_info_s.data_srt)';
    der_fields_t = fieldnames(sort_info_t.data_srt)';

    N_src = size(src_info.r(:,:), 2);
    N_targ = size(targ_info.r(:,:), 2);

    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.

    K_grid2grid = kern_0(grid_info, grid_info);
    K_grid2grid(1:1+size(K_grid2grid,1):end) = 0;
    % Sort the cols of A_spread_s to match the sorted source points
    A_spread_s = A_spread_s(:, sort_info_s.ptid_srt);
    % Sort the cols of A_spread_t to match the sorted target points
    A_spread_t = A_spread_t(:, sort_info_t.ptid_srt);

    AKA = A_spread_t.' * K_grid2grid * A_spread_s;

    % TODO: correct formula for number of corrections
    ncor = grid_info.n_nbr*mean([N_targ,N_src]);

    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.
    iid = zeros(1,ncor);
    jid = zeros(1,ncor);
    vals = zeros(1,ncor);
    id_start = 0;


    % Loop through all of the bins
    for i = 1:size(sort_info_s.id_start, 2) -1
        bin_idx = i - 1; % Because bins are 0-indexed

        % Target points in bin i
        idx_ti_start = sort_info_t.id_start(i);
        idx_ti_end = sort_info_t.id_start(i + 1) - 1;
        
        % targ_pts_in_i = sort_info_t.r_srt(:, idx_ti_start:idx_ti_end);
        targ_info_in_i = [];
        for field = der_fields_t
            targ_info_in_i.(field{1}) = sort_info_t.data_srt.(field{1})(:,idx_ti_start:idx_ti_end);
        end

        % Find the intersecting bins
        intersecting_bin_idxes = intersecting_bins_2d(bin_idx, grid_info, proxy_info);


        for j = 1:size(intersecting_bin_idxes, 2)
            % This iter of the loop does interaction between target bin i and 
            % source bin j
            source_bin_idx = intersecting_bin_idxes(j);

            % Source points in bin j
            idx_sj_start = sort_info_s.id_start(source_bin_idx + 1);
            idx_sj_end = sort_info_s.id_start(source_bin_idx + 2) - 1;

            % src_pts_in_j = sort_info_s.r_srt(:, idx_sj_start:idx_sj_end);
            src_pts_in_j = [];
            for field = der_fields_s
                src_pts_in_j.(field{1}) = sort_info_s.data_srt.(field{1})(:,idx_sj_start:idx_sj_end);
            end

            % Exact near-field interactions
            % K_src_to_targ = kern_st(struct('r', src_pts_in_j), ...
            %                     struct('r', targ_pts_in_i));
            K_src_to_targ = kern_st(src_pts_in_j, ...
                                targ_info_in_i);

            % Compute the approx near-field interactions that should be 
            % subtracted
            AKA_chunk = AKA(idx_ti_start:idx_ti_end, idx_sj_start:idx_sj_end);

            Aloc = K_src_to_targ - AKA_chunk;

            % Arrange local corrections intro a sparse matrix
            is = (idx_ti_start:idx_ti_end);
            js = idx_sj_start:idx_sj_end;
            is = repmat(is(:), 1, size(js,2));
            js = repmat(js(:).', size(is,1), 1);
            n_sparse = numel(Aloc);

            iid(id_start + (1:n_sparse)) = is(:).';
            jid(id_start + (1:n_sparse)) = js(:).';
            vals(id_start + (1:n_sparse)) = Aloc(:).';
            id_start = id_start + n_sparse;
        end

    end
    A_addsub = sparse(iid(1:id_start), jid(1:id_start), vals(1:id_start), N_targ, N_src);

    % Reorder the rows to match the original target point ordering
    A_addsub(sort_info_t.ptid_srt, :) = A_addsub;

    % % Reorder the columns to match the original source point ordering
    A_addsub(:, sort_info_s.ptid_srt) = A_addsub;

end
