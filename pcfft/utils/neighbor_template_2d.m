function [pts, bin_idx] = neighbor_template_2d(grid_info, proxy_info)
    % First, find how many bins are intersecting.
    [int_idx_x, int_idx_y] = intersecting_bins_2d(0, grid_info, proxy_info);
    % Now find a bin s.t. all of the intersecting bins have idx >= 0.
    offset_x = ceil((length(int_idx_x) - 1) / 2);
    offset_y = ceil((length(int_idx_y) - 1) / 2);
    bin_idx = offset_x * grid_info.nbin(2) + offset_y;
    [~, ~, intersecting_bin_idxes] = intersecting_bins_2d(bin_idx, grid_info, proxy_info);

    % For each intersecting bin, get the grid points
    nbins = length(intersecting_bin_idxes);
    npts = nbins * grid_info.nspread^2;
    npts_per_bin = grid_info.nspread^2;
    pts = zeros(2, npts) ;
    for i = 1:length(intersecting_bin_idxes)
        b = intersecting_bin_idxes(i);
        [bin_pts, ~, ~] = grid_pts_for_box_2d(b, grid_info);
        pts(:, (i-1)*npts_per_bin + 1:i*npts_per_bin) = bin_pts;
    end

    % Center at center of bin bin_idx
    % ctr = bin_center(bin_idx, grid_info);
    % pts = pts - ctr;
end