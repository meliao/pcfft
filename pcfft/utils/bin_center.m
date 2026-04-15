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

    if grid_info.dim == 2
        id_y = mod(bin_idx, N_y_bins);
        id_x = floor((bin_idx - id_y) / N_y_bins);
    else
        N_z_bins = grid_info.nbin(3);
        id_z = mod(bin_idx, N_z_bins);
        id_y = mod(floor(bin_idx / N_z_bins), N_y_bins);
        id_x = floor(bin_idx / (N_y_bins * N_z_bins));
    end

    % Compute the first and last gridpoints and find their average

    x_center = Lbd(1) - offset + id_x * dx * nbinpts + (nspread-1)/2 * dx;
    y_center = Lbd(2) - offset + id_y * dx * nbinpts + (nspread-1)/2 * dx;

    if grid_info.dim == 2
        center = [x_center; y_center];
    else
        z_center = Lbd(3) - offset + id_z * dx * nbinpts + (nspread-1)/2 * dx;

        center = [x_center; y_center; z_center];
    end
end
