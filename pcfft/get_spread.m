function [A_spread, sort_info, K_src_to_reg] = get_spread(kern_0, kern_der, ...
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


    % disp("get_spread: id_start:")
    % disp(id_start)

    % We only need to compute the K_reg_to_proxy once, so we
    % will do it here.
    if dim == 2

        [pts_0, center_0, row_idxes_0] = grid_pts_for_box_2d(0, grid_info);
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
    r_local = zeros(dim, size(src_info.r(:,:), 2));
    for i = 1:size(id_start, 2) -1
        

        idx_start = id_start(i);
        idx_end = id_start(i+1) - 1;

        % Get the regular grid points and centers of bin i
        % if dim == 2
        %     [pts_i, center_i, row_idxes_i] = grid_pts_for_box_2d(i-1, grid_info);
        % else
        %     [pts_i, center_i, row_idxes_i] = grid_pts_for_box_3d(i-1, grid_info);
        % end
        center_i = bin_center(i-1, grid_info);
        % disp("get_spread: In bin " + int2str(i))
        % disp("get_spread: idx_start " + int2str(idx_start))
        % disp("get_spread: idx_end " + int2str(idx_end))
        % disp("get_spread: center_i:")
        % disp(center_i)
        % disp("get_spread: grid pts in bin:")
        % disp(pts_i)
        src_pts_in_i = r_sorted(:, idx_start:idx_end);
        % disp("get_spread: src_pts_in_i")
        src_pts_in_i_centered = src_pts_in_i - center_i;

        r_local(:, idx_start:idx_end) = src_pts_in_i_centered;
    end
    src_local = sort_info.data_srt;
    src_local.r = r_local;

    % Compute one whole big K_src_to_proxy, and later we'll 
    % index its rows. K_src_to_proxy has shape (n_proxy, n_src)
    K_src_to_proxy = kern_der_pxy(src_local, proxy_info);
    % K_src_to_reg = K_reg_to_proxy \ K_src_to_proxy;

    % Solve the problem via pseudoinverse
    K_reg_to_proxy_pinv = pinv(K_reg_to_proxy, proxy_info.tol);
    K_src_to_reg = K_reg_to_proxy_pinv * K_src_to_proxy;

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
    

    % Now, loop through the bins and start to fill in A
    % Remember, we 0-indexed the bin IDs
    for i = 0:size(id_start,2) - 2
        idx_start = opdim*(id_start(i+1)-1) + 1;
        idx_end = opdim*(id_start(i+2)-1);
        if idx_end<idx_start, continue, end

        if dim == 2
            [~, ~, row_idxes_i] = grid_pts_for_box_2d(i, grid_info);
        else
            row_idxes_i = grid_ids_for_box_3d(i, grid_info);
        end



        % Do the logging if there is a nonempty set of src point indices
        % if idx_end >= idx_start
        %     disp("get_spread: bin i: " + int2str(i))
        %     disp("get_spread: row_idxes_i:")
        %     disp(row_idxes_i)
        %     disp("get_spread: n_grid_pts: " + int2str(n_grid_pts))
        % else
        %     disp("get_spread: bin i : " + int2str(i) + " empty")
        %     disp("get_spread: idx_start: " +  num2str(idx_start))
        %     disp("get_spread: idx_end: " +  num2str(idx_end))
        % end


        % K_src_to_proxy_i = K_src_to_proxy(:, idx_start:idx_end);

        % block_content = K_reg_to_proxy_pinv * K_src_to_proxy_i;
        % block_content = K_reg_to_proxy \ K_src_to_proxy_i;
        block_content = K_src_to_reg(:,idx_start:idx_end);

        % A_spread(row_idxes_i, idx_start:idx_end) = A_spread(row_idxes_i, idx_start:idx_end) + block_content;
        % A_spread(row_idxes_i, idx_start:idx_end) = block_content;

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

    A_spread = sparse(iid, jid, vals, n_grid_pts, opdim*size(src_info.r(:,:), 2));

    sorted_idxes = opdim*(sorted_idxes-1) + (1:opdim).';
    % Undo the sorting
    A_spread(:, sorted_idxes(:)) = A_spread;

end