function [nbr_binids, nbr_gridpts, nbr_grididxes] = neighbor_template_2d(grid_info, proxy_info, bin_idx)
  % neighbor_template_2d  Neighbor bin template for 2D spreading.
  %
  % Computes the set of neighboring bins and their spreading grid points for a
  % given source bin by shifting the abstract spreading template built for
  % grid_info.center_bin.
  %
  % Inputs
  % ------
  % grid_info  : GridInfo
  %   Regular grid and bin structure.
  % proxy_info : ProxyInfo (or equivalent)
  %   Proxy geometry used to determine which bins intersect.
  % bin_idx    : int
  %   Linear index of the source bin.
  %
  % Outputs
  % -------
  % nbr_binids    : array [1, n_nbr]
  %   Linear bin indices of the neighboring bins. -1 entries mark invalid bins.
  % nbr_gridpts   : array [2, n_nbr_pts]
  %   Coordinates of spreading grid points across all neighbor bins.
  %   Out-of-bounds points are retained with a dummy index.
  % nbr_grididxes : array [1, n_nbr_pts]
  %   Linear grid indices for each column of nbr_gridpts. Out-of-bounds points
  %   are assigned dummy_idx = ngrid(1)*ngrid(2) + 1.

    % Get the abstract spreading template built for center_bin.
    [~, tmpl_pts, tmpl_idxes] = abstract_neighbor_spreading_2D(grid_info, proxy_info);

    % Get neighboring bin IDs (needed by callers such as get_addsub).
    [~, ~, nbr_binids] = intersecting_bins_2d(bin_idx, grid_info, proxy_info);

    % Compute the grid-index shift from center_bin to bin_idx.
    cx = floor(grid_info.center_bin / grid_info.nbin(2));
    cy = mod(grid_info.center_bin, grid_info.nbin(2));
    bx = floor(bin_idx / grid_info.nbin(2));
    by = mod(bin_idx, grid_info.nbin(2));
    delta_ix = (bx - cx) * grid_info.nbinpts;
    delta_iy = (by - cy) * grid_info.nbinpts;

    % Shift template indices to the target bin.
    shifted_idxes = tmpl_idxes + [delta_ix; delta_iy];

    % Determine in-bounds mask and compute linear grid indices.
    ngrid = grid_info.ngrid;
    in_bounds = shifted_idxes(1,:) >= 1 & shifted_idxes(1,:) <= ngrid(1) & ...
                shifted_idxes(2,:) >= 1 & shifted_idxes(2,:) <= ngrid(2);
    dummy_idx = ngrid(1) * ngrid(2) + 1;
    nbr_grididxes = (shifted_idxes(1,:) - 1) * ngrid(2) + shifted_idxes(2,:);
    nbr_grididxes(~in_bounds) = dummy_idx;

    % Compute physical coordinates: tmpl_pts are offsets from center_bin's box
    % center, so adding bin_idx's box center gives absolute coordinates.
    [~, bin_center] = grid_pts_for_box_2d(bin_idx, grid_info);
    nbr_gridpts = tmpl_pts + bin_center;
end
