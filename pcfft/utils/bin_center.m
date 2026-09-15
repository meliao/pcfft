function center = bin_center(bin_idx, grid_info)
% Returns the center(s) of the spreading bin(s) at index <bin_idx>
%
%   The bin indexing follows the scheme used in grid_pts_for_box_2d:
%   bin_idx = id_x * N_y_bins + id_y
%
%   <bin_idx> may be a scalar or an array; the result is [dim, numel(bin_idx)].

    if nargin < 2
        error('bin_center requires bin_idx and grid_info');
    end

    % Hoist property reads: classdef access is slow, so never call this per bin.
    dx = grid_info.dx;
    Lbd = grid_info.Lbd;
    nbin = grid_info.nbin;
    nspread = grid_info.nspread;
    nbinpts = grid_info.nbinpts;
    offset = grid_info.offset;
    dim = grid_info.dim;

    b = double(bin_idx(:).');
    step = dx * nbinpts;
    half = (nspread - 1) / 2 * dx;

    N_y_bins = nbin(2);

    if dim == 2
        id_y = mod(b, N_y_bins);
        id_x = (b - id_y) / N_y_bins;

        center = [Lbd(1) - offset + id_x * step + half; ...
                  Lbd(2) - offset + id_y * step + half];
    else
        N_z_bins = nbin(3);
        id_z = mod(b, N_z_bins);
        id_y = mod(floor(b / N_z_bins), N_y_bins);
        id_x = floor(b / (N_y_bins * N_z_bins));

        center = [Lbd(1) - offset + id_x * step + half; ...
                  Lbd(2) - offset + id_y * step + half; ...
                  Lbd(3) - offset + id_z * step + half];
    end
end
