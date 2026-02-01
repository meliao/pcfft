function [A_spread, K_src_to_reg, sort_info] = get_spread(kern_0, kern, ...
                                            src_info, grid_info, proxy_info)
    % This routine returns the matrix that maps charge strengths at srcinfo.r to 
    % charge strengths on the equispaced grid.
    % Inputs:
    %   kern_0: function handle with calling sequence kern_0(src,targ) 
    %           this is the free-space kernel
    %   kern: function handle with calling sequence kern(src,targ) 
    %           this kernel can have derivatives on the target
    %   src_info: struct with field r (dim, nsrc) source points
    %   grid_info: GridInfo object describing the regular grid
    %   proxy_info: ProxyInfo object describing the proxy points
    % Outputs:
    %   A_spread: sparse matrix of shape (ngrid^dim, nsrc) mapping source
    %             strengths to grid strengths
    %   K_src_to_reg: matrix of shape (nreg, nsrc). Entry (i, j) is the kernel
    %                 evaluation kern_0(z_i - x_j) where z_i is a regular grid 
    %                 point and x_j is a source point.
    %   sort_info: SortInfo object describing the sorting of source points into bins
    dim = proxy_info.dim;


    % First, sort the points into bins
    if dim == 2

        sort_info = SortInfo(src_info.r, grid_info.dx, grid_info.Lbd, ...
                            grid_info.nbin, grid_info.nbinpts);
        r_sorted = sort_info.r_srt;
        sorted_idxes = sort_info.ptid_srt;
        id_start = sort_info.id_start;
    else
        [r_sorted, bin_idxes, id_start] = bin_pts_3d();

    end


    % disp("get_spread: id_start:")
    % disp(id_start)

    % We only need to compute the K_reg_to_proxy once, so we
    % will do it here.
    if dim == 2

        [pts_0, center_0, row_idxes_0] = grid_pts_for_bin_2d(0, grid_info);
    else
        [pts_0, center_0] = grid_pts_for_bin_3d();
    end
    pts_0_centered = pts_0 - center_0;
    K_reg_to_proxy = kern_0(struct('r',pts_0_centered), proxy_info);
    % K_reg_to_proxy_pinv = pinv(K_reg_to_proxy);

    % A is a sparse matrix with shape (ngrid^2, nsrc)
    n_grid_pts = size(grid_info.r, 2);
    A_spread = sparse(n_grid_pts, size(src_info.r, 2));
    % disp("get_spread: A_spread shape: ")
    % disp(size(A_spread))

    % First, loop through the bins and construct "local" source
    % points which are (src points) - (bin center)

    % disp("get_spread: id_start")
    % disp(id_start)
    r_local = zeros(dim, size(src_info.r, 2));
    for i = 1:size(id_start, 2) -1
        

        idx_start = id_start(i);
        idx_end = id_start(i+1) - 1;

        % Get the regular grid points and centers of bin i
        if dim == 2
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_2d(i-1, grid_info);
        else
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_3d();
        end
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


    % Compute one whole big K_src_to_proxy, and later we'll 
    % index its rows. K_src_to_proxy has shape (n_proxy, n_src)
    K_src_to_proxy = kern(struct('r',r_local), proxy_info);
    K_src_to_reg = K_reg_to_proxy \ K_src_to_proxy;

    % Now, loop through the bins and start to fill in A
    % Remember, we 0-indexed the bin IDs
    for i = 0:size(id_start,2) - 2
        if dim == 2
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_2d(i, grid_info);
        else
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_3d();
        end

        idx_start = id_start(i+1);
        idx_end = id_start(i+2) - 1;

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

        A_spread(row_idxes_i, idx_start:idx_end) = A_spread(row_idxes_i, idx_start:idx_end) + block_content;

    end

    % Undo the sorting
    A_spread(:, sorted_idxes) = A_spread;

end