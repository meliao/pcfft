function [id_xs, id_ys, binids] = intersecting_bins_2d(bin_idx, grid_info)
    % Given a set of bins which are described by grid_info, and a set of proxy
    % surfaces which are described by proxy_info, return id_xs and id_ys. The
    % product of these two sets of bins is the set of intersecting bin idxes.
    %
    % This function may return invalid bin idxes in the first two return values,
    % i.e. < 0 or >= grid_info.nbin(d). In the third return value, these invalid
    % bin idxes are set to -1.
    % We say that two bins intersect if their proxy surfaces intersect at all.
    %
    % <bin_idx> may be scalar (outputs are rows, as before) or an array (outputs are [n_nbr, nb]).

    % Hoist property reads: classdef access is slow, so pass all bin indices at once.
    nbin = grid_info.nbin;
    offsets = grid_info.nbr_offsets;

    N_x_bins = nbin(1);
    N_y_bins = nbin(2);

    b = double(bin_idx(:).');
    id_y = mod(b, N_y_bins);
    id_x = (b - id_y) / N_y_bins;

    id_xs = offsets(1, :).' + id_x;   % [n_nbr, nb]
    id_ys = offsets(2, :).' + id_y;

    binids = id_ys + id_xs * N_y_bins;
    binids(id_xs < 0 | id_xs >= N_x_bins | ...
           id_ys < 0 | id_ys >= N_y_bins) = -1;

    if isscalar(bin_idx)
        id_xs = id_xs.';
        id_ys = id_ys.';
        binids = binids.';
    end
end
