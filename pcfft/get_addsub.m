function [A_addsub] = get_addsub(kern_0, kern_st, src_info, ...
    targ_info, grid_info, proxy_info, sort_info_s, sort_info_t, ...
    A_spread_s, A_spread_t)
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
    % src_info : point_info
    %   Specifies the source points.
    % targ_info : point_info
    %   Specifies the target points.
    % grid_info : GridInfo
    %   GridInfo object describing the regular grid.
    % proxy_info : ProxyInfo
    %   ProxyInfo object describing the proxy points.
    % sort_info_s : SortInfo
    %   Specifies how source points are sorted into bins.
    % sort_info_t : SortInfo
    %   Specifies how target points are sorted into bins.
    % A_spread_s : sparse matrix [nreg, opdim(2)*nsrc]
    %   Maps source strengths to equivalent strengths on the regular grid.
    % A_spread_t : sparse matrix [nreg, opdim(1)*ntarg]
    %   Maps target strengths to equivalent strengths on the regular grid. Its adjoint is used to map regular grid strengths to equivalent strengths on the target points.
    %
    % Returns
    % -------
    % A_addsub : sparse matrix [n_targ, n_src]
    %   A sparse matrix which corrects for the incorrect near-field interactions computed using the spreading matrices.

    der_fields_s = fieldnames(sort_info_s.data_srt)';
    der_fields_t = fieldnames(sort_info_t.data_srt)';

    N_src = size(src_info.r(:,:), 2);
    N_targ = size(targ_info.r(:,:), 2);

    dim = size(src_info.r, 1);
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
    if dim == 2
        % [nbr_binids, reg_neighbor_template_pts, ~, nbr_bin_idx] = neighbor_template_2d(grid_info, proxy_info);
        % [pts0, ctr_0, ~] = grid_pts_for_box_2d(nbr_bin_idx, grid_info);
        [pts0, reg_neighbor_template_pts] = abstract_neighbor_spreading_2D(grid_info, proxy_info);
    else
        [~, reg_neighbor_template_pts, ~, nbr_bin_idx] = neighbor_template_3d(grid_info, proxy_info);
        [pts0, ctr_0, ~] = grid_pts_for_box_3d(nbr_bin_idx, grid_info);
    end
    % disp("get_addsub: nbr_binids size: " + int2str(size(nbr_binids)));
    % disp("get_addsub: nbr_binids: ");
    % disp(nbr_binids);
    % pts0_centered = pts0 - ctr_0;
    nbr_info = struct('r', reg_neighbor_template_pts);

    bin_info = struct('r', pts0);
    K_nbr2bin = kern_0(nbr_info, bin_info);
    r = 0;
    for i = 1:dim
        r = r + (reg_neighbor_template_pts(i,:) - pts0(i,:).').^2;
    end
    K_nbr2bin(r<1e-14) = 0;

    % Rows of A_addsub are ordered according to sorted target points.
    % Cols of A_addsub are ordered according to sorted source points.
    % A_addsub = sparse(N_targ, N_src);

    % size of pairwise interaction
    opdim = [size(A_spread_t,2)/N_targ, size(A_spread_s,2)/N_src];

    % Sort the cols of A_spread_s and A_spread_t to match the sorted source points
    src_sort_ids = opdim(2)*(sort_info_s.ptid_srt-1) + (1:opdim(2)).';
    A_spread_s = A_spread_s(:, src_sort_ids(:));
    targ_sort_ids = opdim(1)*(sort_info_t.ptid_srt-1) + (1:opdim(1)).';
    A_spread_t = A_spread_t(:, targ_sort_ids(:));


    % Add 1 row of zeros to A_spread_s to handle empty bins
    A_spread_s = [A_spread_s; sparse(1, opdim(2)*N_src)];
    % dummy_idx = n_gridpts + 1;

    % TODO: correct formula for number of corrections
    ncor = grid_info.n_nbr*ceil(mean([opdim(1)*N_targ,opdim(2)*N_src]));

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
        if dim == 2
            [~, ctr_i, reg_idxs_i] = grid_pts_for_box_2d(bin_idx, grid_info);
        else
            [reg_idxs_i] = grid_ids_for_box_3d(bin_idx, grid_info);
        end

        % Target points in bin i
        idx_ti_start = sort_info_t.id_start(i);
        idx_ti_end = sort_info_t.id_start(i + 1) - 1;
        
        % targ_pts_in_i = sort_info_t.r_srt(:, idx_ti_start:idx_ti_end);
        targ_info_in_i = [];
        for field = der_fields_t
            targ_info_in_i.(field{1}) = sort_info_t.data_srt.(field{1})(:,idx_ti_start:idx_ti_end);
        end

        % Build the spreading template
        if dim == 2
            [nbr_binids, ~, nbr_grididxes, ~] = neighbor_template_2d(grid_info, proxy_info, bin_idx);
            nbr_binids(nbr_binids==-1) =[];
        else
             [nbr_binids, nbr_grididxes] = neighbor_bins_3d(grid_info, proxy_info, bin_idx);
        end

        % Loop through all of the neighbor bins and fill in the local source points. 
        % After this loop, we will update A_add and A_sub with the neigbors of bin i.

        % get index of first and last source in each neighboring bin
        idx_sj_starts = sort_info_s.id_start(nbr_binids + 1);
        idx_sj_ends = sort_info_s.id_start(nbr_binids + 2) - 1;

        % remove empty neighbors
        ifilled = idx_sj_ends>=idx_sj_starts;
        idx_sj_starts = idx_sj_starts(ifilled);
        idx_sj_ends = idx_sj_ends(ifilled);

        % get list of all neighbors
        source_idx = zeros(1,grid_info.n_nbr);
        source_idx_dof = zeros(1,opdim(2)*grid_info.n_nbr);
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
            src_pts_in_j.(field{1}) = sort_info_s.data_srt.(field{1})(:,source_idx);
        end

        % Update A_addsub with exact near-field interactions. This is the "add"
        % part.
        K_src_to_targ = kern_st(src_pts_in_j, ...
                            targ_info_in_i);
        r = 0;
        for k = 1:dim
            r = r + (src_pts_in_j.r(k,:) - targ_info_in_i.r(k,:).').^2;
        end
        r = reshape(r, 1, size(targ_info_in_i.r,2), 1, size(src_pts_in_j.r,2));
        r = repmat(r,opdim(1),1,opdim(2),1);
        r = reshape(r, size(K_src_to_targ));
        K_src_to_targ(r<1e-14) = 0;

        % Update A_sub with approximated near-field interactions. This is the 
        % "sub" part.
        A_spread_t_i = A_spread_t(reg_idxs_i, opdim(1)*(idx_ti_start-1)+1:opdim(1)*idx_ti_end);
        A_spread_s_j = A_spread_s(nbr_grididxes, source_idx_dof);
        % AKA_chunk = (A_spread_t_i.' * K_nbr2bin) * A_spread_s_j;
        disp("get_addsub: size(A_spread_t_i): " + int2str(size(A_spread_t_i)));
        disp("get_addsub: size(K_nbr2bin): " + int2str(size(K_nbr2bin)));
        disp("get_addsub: size(A_spread_s_j): " + int2str(size(A_spread_s_j)));

        AKA_chunk = A_spread_t_i.' * (K_nbr2bin * A_spread_s_j);

        Aloc =  K_src_to_targ - AKA_chunk;

        % Update COO arrays.
        is = (opdim(1)*(idx_ti_start-1)+1:opdim(1)*idx_ti_end);
        js = source_idx_dof;
        is = repmat(is(:), 1, size(js,2));
        js = repmat(js(:).', size(is,1), 1);
        n_sparse = numel(Aloc);

        iid(id_start + (1:n_sparse)) = is(:).';
        jid(id_start + (1:n_sparse)) = js(:).';
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
