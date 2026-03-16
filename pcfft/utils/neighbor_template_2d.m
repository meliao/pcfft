function [nbr_binids, nbr_gridpts, nbr_grididxes, bin_idx] = neighbor_template_2d(grid_info, proxy_info, bin_idx)
  % neighbor_template_2d  Neighbor bin template for 2D spreading.
  %
  % Computes the set of neighboring bins and their spreading grid points for a
  % given source bin. When called without bin_idx, selects an interior reference
  % bin whose full neighborhood lies within the grid.
  %
  % The result serves as a reusable template: for any interior bin, the neighbor
  % structure is the same up to a fixed offset.
  %
  % Inputs
  % ------
  % grid_info  : GridInfo
  %   Regular grid and bin structure.
  % proxy_info : ProxyInfo (or equivalent)
  %   Proxy geometry used to determine which bins intersect.
  % bin_idx    : int, optional
  %   Linear index of the source bin. If omitted, a suitable interior bin is
  %   chosen automatically.
  %
  % Outputs
  % -------
  % nbr_binids    : array [1, n_nbr]
  %   Linear bin indices of the neighboring bins. -1 entries mark invalid bins.
  % nbr_gridpts   : array [2, n_nbr_pts]
  %   Coordinates of unique spreading grid points across all neighbor bins.
  %   Out-of-bounds points are retained with a dummy index.
  % nbr_grididxes : array [1, n_nbr_pts]
  %   Linear grid indices for each column of nbr_gridpts. Out-of-bounds points
  %   are assigned dummy_idx = ngrid(1)*ngrid(2) + 1.
  % bin_idx       : int
  %   The source bin index used.


    if nargin == 2
        % First, find how many bins are intersecting.
        [int_idx_x, int_idx_y] = intersecting_bins_2d(0, grid_info, proxy_info);
        % Filter the invalid ones
        int_idx_x = unique(int_idx_x(int_idx_x >= 0));
        int_idx_y = unique(int_idx_y(int_idx_y >= 0));
        % Now find a bin s.t. all of the intersecting bins have idx >= 0.
        offset_x = length(int_idx_x) - 1;
        offset_y = length(int_idx_y) - 1;
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
    nspread = grid_info.nspread;

    % Preallocate arrays for nbr_gridpts and nbr_grididxes.
    nbr_gridpts = zeros(2, nspread^2 * length(nbr_binids)) * NaN;
    nbr_grididxes = zeros(1, nspread^2 * length(nbr_binids)) * NaN;

    for i = 1:length(nbr_binids)
        % if nbr_binids(i) == -1
        %     continue;
        % end
        [pts, ~, row_idxes] = grid_pts_for_box_2d(nbr_binids(i), grid_info);
        start_idx = (i-1) * nspread^2 + 1;
        end_idx = start_idx + nspread^2 - 1;
        nbr_gridpts(:, start_idx:end_idx) = pts;
        nbr_grididxes(start_idx:end_idx) = row_idxes;
    end

    % Use unique and setdiff to find the unique grid idxes and the corresponding
    % grid pts.
    % TODO
    [~, unique_idx] = unique(nbr_grididxes);
    nbr_grididxes = nbr_grididxes(unique_idx);
    nbr_gridpts = nbr_gridpts(:, unique_idx);

    % Remove all NaNs
    nbr_grididxes = nbr_grididxes(~isnan(nbr_grididxes));
    nbr_gridpts = nbr_gridpts(:, ~isnan(nbr_grididxes));



    % % Build the nbr_gridpts
    % % npts = number of regular grid points across the side of the neighbor template
    % nx = size(int_idx_x, 2);
    % npts = nx * nbinpts + 2 * rpad; 
    % minxnbr = min(int_idx_x); % Minimum x bin index of the neighbor bins
    % minynbr = min(int_idx_y); % Minimum y bin index of the neighbor bins

    % nbr_xpts = Lbd(1) - offset + minxnbr * nbinpts * dx + dx * (0:npts-1);
    % nbr_ypts = Lbd(2) - offset + minynbr * nbinpts * dx + dx * (0:npts-1);
    % [X, Y] = meshgrid(nbr_xpts, nbr_ypts);
    % nbr_gridpts = [X(:).'; Y(:).'];

    % % % Build the nbr_grididxes. This logic is copied from grid_pts_for_box_2d
    % % and npts is used instead of nspread.
    % x_positions =  minxnbr * nbinpts +1 : minxnbr * nbinpts + npts ;
    % y_positions =  minynbr * nbinpts +1  : minynbr * nbinpts + npts ;

    % % Mark out-of-bounds grid points with a dummy index
    ngridpts = grid_info.ngrid(1) * grid_info.ngrid(2);
    dummy_idx =  ngridpts + 1;
    
    % nbr_grididxes = y_positions(:) + (x_positions(:)-1).' * ngrid(2);

    % % Mark the row_idxes corresponding to out-of-bounds grid points with a dummy
    % % Might need a tiny bit of margin here
    margin = 0.1 * dx;

    % Check the y points
    nbr_grididxes(nbr_gridpts(2,:) < rmin(2) - margin | nbr_gridpts(2,:) > rmax(2) + margin,:) = dummy_idx;
    % Check the x points
    nbr_grididxes(:,nbr_gridpts(1,:) < rmin(1) - margin | nbr_gridpts(1,:) > rmax(1) + margin) = dummy_idx;

    % nbr_grididxes = nbr_grididxes(:).';
end