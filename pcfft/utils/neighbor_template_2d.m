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
    rmax = grid_info.rmax;;
    rmin = grid_info.rmin;

    % % Build the nbr_gridpts
    % % npts = number of regular grid points across the side of the neighbor template
    nx = size(int_idx_x, 2);
    npts = nx * nbinpts + 2 * rpad; 
    % % npts_y = length(int_idx_y) * nbinpts + 2 * rpad;
    minxnbr = min(int_idx_x); % Minimum x bin index of the neighbor bins
    minynbr = min(int_idx_y); % Minimum y bin index of the neighbor bins

    nbr_xpts = Lbd(1) - offset + minxnbr * dx + dx * (0:npts-1);
    nbr_ypts = Lbd(2) - offset + minynbr * dx + dx * (0:npts-1);
    [X, Y] = meshgrid(nbr_xpts, nbr_ypts);
    nbr_gridpts = [X(:).'; Y(:).'];

    % % Build the nbr_grididxes. This logic is copied from grid_pts_for_box_2d
    % and npts is used instead of nspread.
    x_positions = minxnbr * nbinpts + 1 : minxnbr * nbinpts + npts;
    y_positions = minynbr * nbinpts + 1 : minynbr * nbinpts + npts;


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


    % Some of the nbr_grididxes may be negative. Flip these.
    % neg_idxes = nbr_grididxes < 1;
    % nbr_grididxes(neg_idxes) = -1 * nbr_grididxes(neg_idxes);

    % Loop through nbr_binids to build nbr_gridpts and nbr_grididxes
    % nbr_gridpts = [];
    % nbr_grididxes = [];

    % nboxpts = (grid_info.nspread)^2;
    % boxidxes = 1:nboxpts;

    % for i = 1:length(nbr_binids)
    %     b = nbr_binids(i);
    %     disp("neighbor_template_2d: Processing neighbor bin id " + int2str(b));
    %     % Get grid pts and grid idxes for box b
    %     [box_pts, ~, box_grididxes] = grid_pts_for_box_2d(b, grid_info);

    %     % disp("neighbor_template_2d: box_grididxes: ")
    %     % disp(box_grididxes);
    %     % Which grididxes already exist in nbr_grididxes?
    %     [~, ia, ~] = intersect(nbr_grididxes, box_grididxes);
    %     % disp("neighbor_template_2d: ia: ")
    %     % disp(ia);

    %     % Append only the new grid pts and grid idxes
    %     new_idxes = setdiff(boxidxes, ia);
    %     % disp("neighbor_template_2d: new_idxes: ")
    %     % disp(new_idxes);
    %     nbr_gridpts = [nbr_gridpts, box_pts(:, new_idxes)];
    %     nbr_grididxes = [nbr_grididxes, box_grididxes(new_idxes)];
    % end


    % Center at center of bin bin_idx
    % ctr = bin_center(bin_idx, grid_info);
    % pts = pts - ctr;
end