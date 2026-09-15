function [A_addsub] = get_addsub(kern_0, kern_st, grid_info, proxy_info, ...
    sort_info_s, sort_info_t, spread_blk_s, spread_blk_t)
    % Compute the correction for near-field interactions.
    %
    % Parameters
    % ----------
    % kern_0 : kernel
    %   Free-space kernel, which must be scalar-valued
    % kern_st : kernel
    %   Direct interaction kernel, which must be a combination of
    %   derivatives of the free-space kernel. If left empty the free-space
    %   kernel will be used.
    %   Each pairwise interaction must of shape [opdim(1), opdim(2)]
    % grid_info : GridInfo
    %   GridInfo object describing the regular grid.
    % proxy_info : ProxyInfo
    %   ProxyInfo object describing the proxy points.
    % sort_info_s : SortInfo
    %   Specifies how source points are sorted into bins.
    % sort_info_t : SortInfo
    %   Specifies how target points are sorted into bins.
    % spread_blk_s : matrix [nspread^dim, opdim(1)*nsrc]
    %   Dense source spreading weights in sorted-point order, the third output
    %   of ``get_spread()`` for the sources.
    % spread_blk_t : matrix [nspread^dim, opdim(1)*ntarg]
    %   Dense target spreading weights in sorted-point order, the third output
    %   of ``get_spread()`` for the targets.
    %
    % Returns
    % -------
    % A_addsub : sparse matrix [n_targ, n_src]
    %   A sparse matrix which corrects for the incorrect near-field interactions computed using the spreading matrices.

    der_fields_s = fieldnames(sort_info_s.data_srt)';
    der_fields_t = fieldnames(sort_info_t.data_srt)';

    N_src = size(sort_info_s.r_srt(:,:), 2);
    N_targ = size(sort_info_t.r_srt(:,:), 2);

    dim = size(sort_info_s.r_srt, 1);
    if isempty(kern_st), kern_st = kern_0; end

    if ~isa(kern_0,'function_handle')
        try
            kern_0 = kern_0.eval;
        catch
            error('kern_0 is not a function and does not have an eval property')
        end
    end
    if ~isa(kern_st,'function_handle')
        try
            kern_st = kern_st.eval;
        catch
            error('kern_st is not a function and does not have an eval property')
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build a spreading template matrix for adjacent source points.
    % Then build a list of regular gridpoints that are in the intersecting bins
    % template_pos(:, k) gives the columns of the spreading template covered by
    % the spreading box of the k-th entry of grid_info.nbr_offsets, in the same
    % point order as the rows of spread_blk_s / spread_blk_t.
    if dim == 2
        [pts0, reg_neighbor_template_pts, ~, template_pos] = abstract_neighbor_spreading_2D(grid_info, proxy_info);
    else
        [pts0, reg_neighbor_template_pts, ~, template_pos] = abstract_neighbor_spreading_3D(grid_info, proxy_info);
    end
    box_center = bin_center(grid_info.center_bin, grid_info);
    pts0 = pts0 - box_center;

    nbr_info = struct('r', reg_neighbor_template_pts);

    box_info = struct('r', pts0);
    K_nbr2bin = kern_0(nbr_info, box_info);
    r = 0;
    for i = 1:dim
        r = r + (reg_neighbor_template_pts(i,:) - pts0(i,:).').^2;
    end
    K_nbr2bin(r<1e-14) = 0;

    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.
    % A_addsub = sparse(N_targ, N_src);

    % size of pairwise interaction
    opdim = [size(spread_blk_t,2)/N_targ, size(spread_blk_s,2)/N_src];

    src_sort_ids = opdim(2)*(sort_info_s.ptid_srt-1) + (1:opdim(2)).';
    targ_sort_ids = opdim(1)*(sort_info_t.ptid_srt-1) + (1:opdim(1)).';

    id_start = 0;

    nbins_total = size(sort_info_s.id_start, 2) - 1;
    nbin_l = grid_info.nbin;
    offsets_l = grid_info.nbr_offsets;
    bb_all = 0:nbins_total-1;
    if dim == 2
        Ny = nbin_l(2);
        bin_idy = mod(bb_all, Ny);
        bin_idx_x = (bb_all - bin_idy) / Ny;
        lin_off = offsets_l(1, :) * Ny + offsets_l(2, :);
        off_x = offsets_l(1, :);
        off_y = offsets_l(2, :);
    else
        Ny = nbin_l(2);
        Nz = nbin_l(3);
        bin_idz = mod(bb_all, Nz);
        bin_idy = mod(floor(bb_all / Nz), Ny);
        bin_idx_x = floor(bb_all / (Ny * Nz));
        lin_off = offsets_l(1, :) * (Ny * Nz) + offsets_l(2, :) * Nz + offsets_l(3, :);
        off_x = offsets_l(1, :);
        off_y = offsets_l(2, :);
        off_z = offsets_l(3, :);
    end
    Nx = nbin_l(1);

    % Fetch the remaining per-bin properties
    s_id_start = sort_info_s.id_start;
    t_id_start = sort_info_t.id_start;
    s_data_srt = sort_info_s.data_srt;
    t_data_srt = sort_info_t.data_srt;
    n_nbr = grid_info.n_nbr;

    % "sub" strategy: template contracts vs the whole stencil, direct per-neighbour (better when bins are mostly empty).
    nbox_l = size(K_nbr2bin, 1);
    ntemplate_l = size(K_nbr2bin, 2);

    % Count occupied neighbour bins per occupied target bin (one pass per offset).
    s_counts = diff(s_id_start);
    t_counts = diff(t_id_start);
    n_occ_nbr = 0;
    npair = 0;
    for kk = 1:size(offsets_l, 2)
        okk = off_x(kk) + bin_idx_x >= 0 & off_x(kk) + bin_idx_x < Nx & ...
              off_y(kk) + bin_idy >= 0 & off_y(kk) + bin_idy < Ny;
        if dim == 3
            okk = okk & off_z(kk) + bin_idz >= 0 & off_z(kk) + bin_idz < Nz;
        end
        idxk = find(okk);
        ct = double(t_counts(idxk));
        cs = double(s_counts(idxk + lin_off(kk)));
        n_occ_nbr = n_occ_nbr + nnz(ct > 0 & cs > 0);
        npair = npair + sum(ct .* cs);
    end
    mean_occ_nbr = n_occ_nbr / max(1, nnz(t_counts > 0));

    use_template = (nbox_l * mean_occ_nbr >= ntemplate_l);

    ncor = npair * opdim(1) * opdim(2);

    % These are the arrays we will use to build the sparse A_addsub
    % in COO format.
    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.
    iid = zeros(1,ncor);
    jid = zeros(1,ncor);
    vals = zeros(1,ncor);

    % Template path: C is bin-independent, so compute it for a block of bins at once.
    max_block_bytes = 2^28;   % 256 MB
    bytes_per_entry = 16;     % complex double
    max_block_rows = max(1, floor(max_block_bytes / ...
                         (size(K_nbr2bin,2) * bytes_per_entry)));
    block_bin_end = 0;    % last bin covered by the current block
    block_row0 = 0;       % row offset of the current block in spread_blk_t cols
    C_block = [];

    % Loop through all of the bins. 
    for i = 1:nbins_total
        bin_idx = i - 1; % Because bins are 0-indexed
        % disp("get_addsub: Processing bin " + int2str(bin_idx));

        % Target points in bin i
        idx_ti_start = t_id_start(i);
        idx_ti_end = t_id_start(i + 1) - 1;
        if idx_ti_start > idx_ti_end, continue, end
        
        % targ_pts_in_i = sort_info_t.r_srt(:, idx_ti_start:idx_ti_end);
        targ_info_in_i = [];
        for field = der_fields_t
            targ_info_in_i.(field{1}) = t_data_srt.(field{1})(:,idx_ti_start:idx_ti_end);
        end

        % Find the neighboring bins of bin i. off_ids records which entry of
        % grid_info.nbr_offsets each surviving neighbor came from, which is
        % what indexes template_pos.
        ix = off_x + bin_idx_x(i);
        iy = off_y + bin_idy(i);
        valid = ix >= 0 & ix < Nx & iy >= 0 & iy < Ny;
        if dim == 3
            iz = off_z + bin_idz(i);
            valid = valid & iz >= 0 & iz < Nz;
        end
        off_ids = find(valid);
        nbr_binids = bin_idx + lin_off(off_ids);

        % Loop through all of the neighbor bins and fill in the local source points. 
        % After this loop, we will update A_add and A_sub with the neigbors of bin i.

        % get index of first and last source in each neighboring bin
        idx_sj_starts = s_id_start(nbr_binids + 1);
        idx_sj_ends = s_id_start(nbr_binids + 2) - 1;

        % remove empty neighbors
        ifilled = idx_sj_ends>=idx_sj_starts;
        idx_sj_starts = idx_sj_starts(ifilled);
        idx_sj_ends = idx_sj_ends(ifilled);
        off_ids = off_ids(ifilled);

        % get list of all neighbors
        source_idx = zeros(1,n_nbr);
        source_idx_dof = zeros(1,opdim(2)*n_nbr);
        istart = 1;
        for j = 1:length(idx_sj_starts)
            % This iter of the loop does interaction between target bin i and 
            % the jth neighboring source bin
            idx_sj_start = idx_sj_starts(j);
            idx_sj_end = idx_sj_ends(j);

            % store indices
            source_idx(istart:istart+(idx_sj_end-idx_sj_start)) = idx_sj_start:idx_sj_end;
            source_idx_dof(opdim(2)*(istart-1)+1:opdim(2)*(istart+idx_sj_end-idx_sj_start)) = opdim(2)*(idx_sj_start-1)+1:opdim(2)*idx_sj_end;
            istart = istart + (idx_sj_end-idx_sj_start+1);
        end
        source_idx = source_idx(1:istart-1);
        source_idx_dof = source_idx_dof(1:opdim(2)*(istart-1));

        % It may be the case that there are no source points in the bins 
        % neighboring target bin i. 
        if isempty(source_idx)
            continue;
        end

        src_pts_in_j = [];
        for field = der_fields_s
            src_pts_in_j.(field{1}) = s_data_srt.(field{1})(:,source_idx);
        end

        % Update A_addsub with exact near-field interactions. This is the "add"
        % part.
        K_src_to_targ = kern_st(src_pts_in_j, ...
                            targ_info_in_i);
        % Zero out the self interactions. 
        r = 0;
        for k = 1:dim
            r = r + (src_pts_in_j.r(k,:) - targ_info_in_i.r(k,:).').^2;
        end
        if opdim(1) == 1 && opdim(2) == 1
            K_src_to_targ(r < 1e-14) = 0;
        else
            [ti, sj] = find(r < 1e-14);
            if ~isempty(ti)
                rows = opdim(1)*(ti(:).'-1) + (1:opdim(1)).';   % [opdim1, npair]
                cols = opdim(2)*(sj(:).'-1) + (1:opdim(2)).';   % [opdim2, npair]
                lin = reshape(rows, opdim(1), 1, []) + ...
                      (reshape(cols, 1, opdim(2), []) - 1) * size(K_src_to_targ,1);
                K_src_to_targ(lin(:)) = 0;
            end
        end

        % Update A_sub with approximated near-field interactions. This is the
        % "sub" part. Both spreading blocks are contiguous column slices of the
        % dense weights returned by get_spread, so no sparse indexing is needed.
        % Refill C_block once this bin runs past it.
        if use_template && i > block_bin_end
            block_row0 = opdim(1)*(t_id_start(i)-1);
            jb = i;
            while jb < nbins_total && ...
                  opdim(1)*(t_id_start(jb+2)-1) - block_row0 <= max_block_rows
                jb = jb + 1;
            end
            block_bin_end = jb;
            block_cols = block_row0+1 : opdim(1)*(t_id_start(block_bin_end+1)-1);
            C_block = spread_blk_t(:, block_cols).' * K_nbr2bin;
        end

        % Row range of this bin; index C_block directly rather than slicing out a copy.
        cols_ti = opdim(1)*(idx_ti_start-1)+1 : opdim(1)*idx_ti_end;

        % Hit each neighbor with just the template columns its box covers.
        if use_template
            rows_i = cols_ti - block_row0;
            AKA_chunk = zeros(numel(rows_i), numel(source_idx_dof), 'like', C_block);
            col = 0;
            for j = 1:length(idx_sj_starts)
                cs = opdim(2)*(idx_sj_starts(j)-1)+1 : opdim(2)*idx_sj_ends(j);
                AKA_chunk(:, col + (1:numel(cs))) = ...
                    C_block(rows_i, template_pos(:, off_ids(j))) * spread_blk_s(:, cs);
                col = col + numel(cs);
            end
        else
            A_spread_t_i = spread_blk_t(:, cols_ti);
            AKA_chunk = zeros(numel(cols_ti), numel(source_idx_dof), 'like', K_nbr2bin);
            col = 0;
            for j = 1:length(idx_sj_starts)
                cs = opdim(2)*(idx_sj_starts(j)-1)+1 : opdim(2)*idx_sj_ends(j);
                AKA_chunk(:, col + (1:numel(cs))) = ...
                    (A_spread_t_i.' * K_nbr2bin(:, template_pos(:, off_ids(j)))) * ...
                    spread_blk_s(:, cs);
                col = col + numel(cs);
            end
        end

        Aloc =  K_src_to_targ - AKA_chunk;

        % Update COO arrays.
        is = opdim(1)*(idx_ti_start-1)+1 : opdim(1)*idx_ti_end;
        js = source_idx_dof;
        nrow = numel(is);
        n_sparse = numel(Aloc);

        if id_start + n_sparse > numel(vals)
            % Extend vectors if necessary
            newlen = max(2*numel(vals), id_start + n_sparse);
            iid(newlen) = 0; jid(newlen) = 0; vals(newlen) = 0;
        end
        iid(id_start + (1:n_sparse)) = repmat(is, 1, numel(js));
        jid(id_start + (1:n_sparse)) = repelem(js, nrow);
        vals(id_start + (1:n_sparse)) = Aloc(:).';
        id_start = id_start + n_sparse;

    end
    iid = iid(1:id_start);
    jid = jid(1:id_start);
    vals = vals(1:id_start);

    targ_sort_ids = targ_sort_ids(:);
    src_sort_ids = src_sort_ids(:);
    iid = targ_sort_ids(iid);
    jid = src_sort_ids(jid);

    % isort = randperm(id_start);
    A_addsub = sparse(iid, jid, vals, opdim(1)*N_targ, opdim(2)*N_src);
    % A_addsub = sparse(iid(isort), jid(isort), vals(isort), opdim(1)*N_targ, opdim(2)*N_src);
    % [jid,isort] = sort(jid);
    % iid = iid(isort);
    % vals = vals(isort);
    % A_addsub = sparse(iid, jid, vals, opdim(1)*N_targ, opdim(2)*N_src);
    % [iid,isort] = sort(iid);
    % jid = jid(isort);
    % vals = vals(isort);
    % 
    % 
    % A_addsub = sparse(iid, jid, vals, opdim(1)*N_targ, opdim(2)*N_src);

    % % Reorder the rows to match the original target point ordering
    % A_addsub(targ_sort_ids, :) = A_addsub;
    % 
    % % % Reorder the columns to match the original source point ordering
    % A_addsub(:, src_sort_ids) = A_addsub; 
end
