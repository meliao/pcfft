function [nbr_binids, nbr_grididxes] = neighbor_bins_3d(grid_info, bin_idx)

    [int_idx_x, int_idx_y, int_idx_z, nbr_binids] = intersecting_bins_3d(bin_idx, grid_info);
    

    % disp("neighbor_template_2d: For bin_idx " + int2str(bin_idx) + ...
    %     ", ind_idx_x: ");
    % disp(ind_idx_x);

    % Info from grid_info that will be used later
    rpad = grid_info.rpad;
    nbinpts = grid_info.nbinpts;
    Lbd = grid_info.Lbd;
    dx = grid_info.dx;
    offset = grid_info.offset;
    ngrid = grid_info.ngrid;
    rmax = grid_info.rmax;
    rmin = grid_info.rmin;

    % % Build the nbr_gridpts
    nx = size(int_idx_x, 2);
    npts = nx * nbinpts + 2 * rpad; 
    minxnbr = min(int_idx_x); % Minimum x bin index of the neighbor bins
    minynbr = min(int_idx_y); % Minimum y bin index of the neighbor bins
    minznbr = min(int_idx_z); % Minimum z bin index of the neighbor bins

    nbr_xpts = Lbd(1) - offset + minxnbr * nbinpts * dx + dx * (0:npts-1);
    nbr_ypts = Lbd(2) - offset + minynbr * nbinpts * dx + dx * (0:npts-1);
    nbr_zpts = Lbd(3) - offset + minznbr * nbinpts * dx + dx * (0:npts-1);

    % % Build the nbr_grididxes. This logic is copied from grid_pts_for_box_3d
    % and npts is used instead of nspread.
    x_positions =  minxnbr * nbinpts +1 : minxnbr * nbinpts + npts ;
    y_positions =  minynbr * nbinpts +1  : minynbr * nbinpts + npts ;
    z_positions =  minznbr * nbinpts +1  : minznbr * nbinpts + npts ;

    % Mark out-of-bounds grid points with a dummy index
    ngridpts = grid_info.ngrid(1) * grid_info.ngrid(2) * grid_info.ngrid(3);
    dummy_idx =  ngridpts + 1;

    nbr_grididxes = (z_positions(:)) + (y_positions(:)-1).'*ngrid(3) + reshape(x_positions-1,1,1,[]) * ngrid(2) * ngrid(3);

    % Mark the row_idxes corresponding to out-of-bounds grid points with a dummy
    % Might need a tiny bit of margin here
    margin = 0.1 * dx;

    nbr_grididxes(nbr_zpts<rmin(3) - margin | nbr_zpts > rmax(3) + margin,:,:) = dummy_idx;
    nbr_grididxes(:,nbr_ypts<rmin(2) - margin | nbr_ypts > rmax(2) + margin,:) = dummy_idx;
    nbr_grididxes(:,:,nbr_xpts<rmin(1) - margin | nbr_xpts > rmax(1) + margin) = dummy_idx;

    nbr_grididxes = nbr_grididxes(:).';
end