function [nbr_binids, nbr_gridpts, nbr_grididxes, bin_idx] = neighbor_template_2d(grid_info, proxy_info, bin_idx)

    if nargin == 2
        % First, find how many bins are intersecting.
        [int_idx_x, int_idx_y] = intersecting_bins_2d(0, grid_info, proxy_info);
        % Now find a bin s.t. all of the intersecting bins have idx >= 0.
        offset_x = ceil((length(int_idx_x) - 1) / 2);
        offset_y = ceil((length(int_idx_y) - 1) / 2);
        bin_idx = offset_x * grid_info.nbin(2) + offset_y;
    end

    [int_idx_x, int_idx_y, nbr_binids] = intersecting_bins_2d(bin_idx, grid_info, proxy_info);
    

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
    % % npts = number of regular grid points across the side of the neighbor template
    nx = size(int_idx_x, 2);
    npts = nx * nbinpts + 2 * rpad; 
    minxnbr = min(int_idx_x); % Minimum x bin index of the neighbor bins
    minynbr = min(int_idx_y); % Minimum y bin index of the neighbor bins

    nbr_xpts = Lbd(1) - offset + minxnbr * nbinpts * dx + dx * (0:npts-1);
    nbr_ypts = Lbd(2) - offset + minynbr * nbinpts * dx + dx * (0:npts-1);
    [X, Y] = meshgrid(nbr_xpts, nbr_ypts);
    nbr_gridpts = [X(:).'; Y(:).'];

    % % Build the nbr_grididxes. This logic is copied from grid_pts_for_box_2d
    % and npts is used instead of nspread.
    x_positions =  minxnbr * nbinpts +1 : minxnbr * nbinpts + npts ;
    y_positions =  minynbr * nbinpts +1  : minynbr * nbinpts + npts ;

    % Mark out-of-bounds grid points with a dummy index
    ngridpts = grid_info.ngrid(1) * grid_info.ngrid(2);
    dummy_idx =  ngridpts + 1;

    nbr_grididxes = ones(1, npts^2);
    for i = 1:npts
        start_row = (i-1) * npts + 1;
        end_row = start_row + npts - 1;

        nbr_grididxes(start_row:end_row) = (x_positions(i)-1) * ngrid(2) + y_positions ;
    end

    % Mark the row_idxes corresponding to out-of-bounds grid points with a dummy
    out_of_bounds = (nbr_gridpts(1, :) < rmin(1)) | (nbr_gridpts(1, :) > rmax(1)) | ...
                    (nbr_gridpts(2, :) < rmin(2)) | (nbr_gridpts(2, :) > rmax(2));
    nbr_grididxes(out_of_bounds) = dummy_idx;


end