function center = bin_center(bin_idx, grid_info)
% Returns the center of the spreading bin at index <bin_idx>
%
%   The bin indexing follows the scheme used in grid_pts_for_box_2d: 
%   bin_idx = id_x * N_y_bins + id_y

    if nargin < 2
        error('bin_center requires bin_idx and grid_info');
    end

    dx = grid_info.dx;
    Lbd = grid_info.Lbd;
    N_y_bins = grid_info.nbin(2);
    nspread = grid_info.nspread;
    nbinpts = grid_info.nbinpts;
    offset = grid_info.offset;

    % Compute integer bin ids (consistent with grid_pts_for_box_2d)
    id_y = mod(bin_idx, N_y_bins);
    id_x = floor((bin_idx - id_y) / N_y_bins);

    % Compute the first and last gridpoints and find their average
    x_first = Lbd(1) - offset + id_x * dx * nbinpts + 0 * dx;
    x_last  = Lbd(1) - offset + id_x * dx * nbinpts + (nspread-1) * dx;

    y_first = Lbd(2) - offset + id_y * dx * nbinpts + 0 * dx;
    y_last  = Lbd(2) - offset + id_y * dx * nbinpts + (nspread-1) * dx;

    center = [(x_first + x_last) / 2; (y_first + y_last) / 2];
end
