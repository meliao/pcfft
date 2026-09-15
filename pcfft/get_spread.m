function [A_spread, sort_info, spread_blk] = get_spread(kern_0, kern_der, ...
                                            src_info, grid_info, proxy_info, der_fields)
    % This routine returns the matrix that maps charge strengths at srcinfo.r to 
    % charge strengths on the equispaced grid.
    %
    % Parameters
    % ----------
    % kern_0 : kernel
    %   The free-space kernel, which must be scalar-valued
    % kern_der : kernel
    %   Some derivative of the free-space kernel. This can be left empty by 
    %   passing in an empty array, in which case the free-space kernel will be used.
    %   Each pairwise interaction must of shape [1, opdim]
    % src_info : point_info
    %   Specifies the source points.
    % grid_info : GridInfo
    %   object describing the regular grid
    % proxy_info : ProxyInfo
    %   object describing the proxy points
    % der_fields : cell array, optional
    %   Contains field names that must be attached to the source point in
    %   ``kern_der``, e.g. {'r', 'n', 'kappa'}. If omitted or left
    %   empty this argument defaults to {'r'}.
    %
    %
    % Returns
    % -------
    % A_spread : sparse matrix [nreg, opdim*nsrc]
    %   Maps source strengths to equivalent strengths on the regular grid.
    % sort_info : SortInfo
    %   Object describing the sorting of source points into bins
    % spread_blk : matrix [nspread^dim, opdim*nsrc]
    %   Dense spreading weights, before they are scattered into A_spread.
    %   This matrix should be passed to get_addsub. Column blocks are ordered 
    %   by *sorted* point index. For more information see the developer note 
    %   in the source code.


    % DEVELOPER NOTE:
    % Suppose point with sorted index p lives in bin i. Bin i has some 
    % associated spreading box. Row p of spread_blk contains the spreading 
    % weights to map point p onto the points of bin i's spreading box.
    % This is meant to make it easy for get_addsub to use the spreading 
    % weights.

    if nargin < 6; der_fields = {}; end
    dim = proxy_info.dim;

    % if kern_der is not provided, we use the free-space kernel
    if isempty(kern_der), kern_der = kern_0; end

    if ~isa(kern_0,'function_handle')
        try
            kern_0 = kern_0.eval;
        catch
            error('kern_0 is not a function and does not have an eval property')
        end
    end
    if ~isa(kern_der,'function_handle')
        try
            kern_der = kern_der.eval;
        catch
            error('kern_der is not a function and does not have an eval property')
        end
    end

    proxy_der = proxy_info.proxy_der;
    kern_0_pxy = @(s,t) wrap_kern_der(kern_0, s, t, proxy_der);
    kern_der_pxy = @(s,t) wrap_kern_der(kern_der, s, t, proxy_der);
    % First, sort the points into bins
    sort_info = SortInfo(src_info, grid_info.dx, grid_info.Lbd, ...
                        grid_info.nbin, grid_info.nbinpts,der_fields);
    r_sorted = sort_info.r_srt;
    sorted_idxes = sort_info.ptid_srt;
    id_start = sort_info.id_start;


    % We only need to compute the K_reg_to_proxy once, so we
    % will do it here.
    if dim == 2

        [pts_0, center_0] = grid_pts_for_box_2d(0, grid_info);
    else
        [pts_0, center_0] = grid_pts_for_box_3d(0, grid_info);
    end
    pts_0_centered = pts_0 - center_0;
    K_reg_to_proxy = kern_0_pxy(struct('r',pts_0_centered), proxy_info);
    if any(size(K_reg_to_proxy) ~= [(proxy_der+1)*size(proxy_info.r,2), size(pts_0_centered,2)])
        error('kern_0 must be scalar-valued. Use kern_component for vector-valued free-space kernels.')
    end
        
    % K_reg_to_proxy_pinv = pinv(K_reg_to_proxy);


    % First, loop through the bins and construct "local" source
    % points which are (src points) - (bin center)

    % disp("get_spread: id_start")
    % disp(id_start)
    % Subtract each source point's bin center.
    nbins = size(id_start, 2) - 1;
    r_local = zeros(dim, size(src_info.r(:,:), 2));
    bin_ctrs = bin_center(0:nbins-1, grid_info);
    npts_per_bin = diff(id_start);
    npts_sorted = id_start(end) - 1;
    if npts_sorted > 0
        ptbin = repelem(1:nbins, npts_per_bin);
        r_local(:, 1:npts_sorted) = r_sorted(:, 1:npts_sorted) - bin_ctrs(:, ptbin);
    end
    src_local = sort_info.data_srt;
    src_local.r = r_local;

    % Compute one whole big K_src_to_proxy, and later we'll 
    % index its rows. K_src_to_proxy has shape (n_proxy, n_src)
    % There are a small number of proxy points so this is not too expensive.
    K_src_to_proxy = kern_der_pxy(src_local, proxy_info);
    % K_src_to_reg = K_reg_to_proxy \ K_src_to_proxy;
    K_src_to_reg = lsqminnorm(K_reg_to_proxy, K_src_to_proxy, proxy_info.tol / 10);
    
    
    % determine dimension of the kernel
    opdim = size(K_src_to_proxy,2)/size(src_info.r(:,:), 2);

    % A is a sparse matrix with shape (ngrid^2, nsrc)
    n_grid_pts = size(grid_info.r, 2);
    % A_spread = sparse(n_grid_pts, opdim*size(src_info.r(:,:), 2));
    % disp("get_spread: A_spread shape: ")
    % disp(size(A_spread))

    num_spread = size(K_src_to_reg,1)*opdim*size(src_info.r(:,:), 2);
    iid = zeros(1,num_spread);
    jid = zeros(1,num_spread);
    vals = zeros(1,num_spread);
    id_id = 0;
    

    % Box row indices are a fixed template plus a per-bin offset; precompute both.
    ngrid = grid_info.ngrid;
    nspread_l = grid_info.nspread;
    nbinpts_l = grid_info.nbinpts;
    nbin_l = grid_info.nbin;
    ii = (1:nspread_l).';
    if dim == 2
        bb = 0:nbins-1;
        id_y_all = mod(bb, nbin_l(2));
        id_x_all = (bb - id_y_all) / nbin_l(2);
        row_base = id_x_all * (nbinpts_l * ngrid(2)) + id_y_all * nbinpts_l;
        row_template = (ii - 1) * ngrid(2) + (1:nspread_l);
        row_template = reshape(row_template.', 1, []);   % i outer, j inner
    else
        bb = 0:nbins-1;
        id_z_all = mod(bb, nbin_l(3));
        id_y_all = mod(floor(bb / nbin_l(3)), nbin_l(2));
        id_x_all = floor(bb / (nbin_l(2) * nbin_l(3)));
        row_base = id_z_all * nbinpts_l + ...
                   id_y_all * (nbinpts_l * ngrid(3)) + ...
                   id_x_all * (nbinpts_l * ngrid(2) * ngrid(3));
        row_template = ii + ((1:nspread_l) - 1) * ngrid(3) + ...
                       reshape((0:nspread_l-1) * ngrid(2) * ngrid(3), 1, 1, []);
        row_template = reshape(row_template, 1, []);
    end

    % Now, loop through the bins and start to fill in A
    % Remember, we 0-indexed the bin IDs
    for i = 0:size(id_start,2) - 2
        idx_start = opdim*(id_start(i+1)-1) + 1;
        idx_end = opdim*(id_start(i+2)-1);
        if idx_end<idx_start, continue, end

        row_idxes_i = row_base(i+1) + row_template;

        block_content = K_src_to_reg(:,idx_start:idx_end);


        % Update COO arrays.
        is = (row_idxes_i);
        js = idx_start:idx_end;
        is = repmat(is(:), 1, size(js,2));
        js = repmat(js(:).', size(is,1), 1);
        n_sparse = numel(block_content);

        iid(id_id + (1:n_sparse)) = is(:).';
        jid(id_id + (1:n_sparse)) = js(:).';
        vals(id_id + (1:n_sparse)) = block_content(:).';
        id_id = id_id + n_sparse;

    end
    iid = iid(1:id_id);
    jid = jid(1:id_id);
    vals = vals(1:id_id);

    % Undo the sorting. Permuting the COO column indices is equivalent to the
    % column assignment A_spread(:, sorted_idxes) = A_spread, but avoids
    % building and then permuting a sparse matrix.
    sorted_idxes = opdim*(sorted_idxes-1) + (1:opdim).';
    sorted_idxes = sorted_idxes(:);
    jid = reshape(sorted_idxes(jid), 1, []);

    A_spread = sparse(iid, jid, vals, n_grid_pts, opdim*size(src_info.r(:,:), 2));

    spread_blk = K_src_to_reg;

end