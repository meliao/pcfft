function A_spread = get_spread(kern_0, kern, src_info, grid_info, proxy_info, nbin)
    % This routine returns the matrix that maps charge strengths at srcinfo.r to 
    % charge strengths on the equispaced grid.
    % We want to 
    if nargin < 6
        nbin = 1;
    end
    dim = proxy_info.dim;


    % First, sort the points into bins
    if dim == 2
        n_grid_pts = (grid_info.ngrid(1) * grid_info.rpad + 1) * (grid_info.ngrid(2) * grid_info.rpad + 1);
        [r_sorted, bin_idxes, id_start] = bin_pts_2d(src_info.r, grid_info.dx, grid_info.ngrid, grid_info.Lbd, nbin);
    else
        [r_sorted, bin_idxes, id_start] = bin_pts_3d(src_info.r, grid_info.dx, grid_info.ngrid, grid_info.Lbd, nbin);
        n_grid_pts = grid_info.ngrid(1) * grid_info.ngrid(2) * grid_info.ngrid(3);

    end

    % We only need to compute the K_reg_to_proxy once, so we
    % will do it here.
    if dim == 2

        [pts_0, center_0, row_idxes_0] = grid_pts_for_bin_2d(0, grid_info, nbin);
    else
        [pts_0, center_0] = grid_pts_for_bin_3d(0, grid_info.dx, grid_info.ngrid, grid_info.Lbd, nbin);
    end
    pts_0_centered = pts_0 - center_0;
    K_reg_to_proxy = kern_0(pts_0_centered, proxy_info.r);
    K_reg_to_proxy_pinv = pinv(K_reg_to_proxy);

    % A is a sparse matrix with shape (ngrid^2, nsrc)
    A_spread = sparse(n_grid_pts, size(src_info.r, 2));

    % First, loop through the bins and construct "local" source
    % points which are (src points) - (bin center)

    disp("id_start")
    disp(id_start)
    r_local = zeros(dim, size(src_info.r, 2));
    for i = 1:size(id_start, 2) -1
        

        idx_start = id_start(i);
        idx_end = id_start(i+1) - 1;

        % Get the regular grid points and centers of bin i
        if dim == 2
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_2d(i, grid_info, nbin);
        else
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_3d();
        end
        disp("In bin " + int2str(i))
        disp("idx_start " + int2str(idx_start))
        disp("idx_end " + int2str(idx_end))
        src_pts_in_i = r_sorted(idx_start:idx_end);
        src_pts_in_i_centered = src_pts_in_i - center_i;

        r_local(:, idx_start:idx_end) = src_pts_in_i_centered;

    end


    % Compute one whole big K_src_to_proxy, and later we'll 
    % index its rows
    K_src_to_proxy = kern_0(r_local, proxy_info.r);


    % Now, loop through the bins and start to fill in A
    % Remember, we 0-indexed the bin IDs
    for i = 1:size(id_start,2)
        if dim == 2
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_2d(i, grid_info, nbin);
        else
            [pts_i, center_i, row_idxes_i] = grid_pts_for_bin_3d();
        end
        sources_in_bin_i = bin_idxes == i;
        K_src_to_proxy_i = K_src_to_proxy(:, sources_in_bin_i);

        block_content = K_reg_to_proxy_pinv * K_src_to_proxy_i;

        A_spread(row_idxes_i, sources_in_bin_i) = block_content;


    end

end