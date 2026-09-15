function [id_xs, id_ys, id_zs, binids] = intersecting_bins_3d(bin_idx, grid_info)
    % Given a set of bins which are described by grid_info, and a set of proxy
    % surfaces which are described by proxy_info, return id_xs, id_ys, id_zs.
    % The product of these sets of bins is the set of intersecting bin idxes.
    %
    % This function may return invalid bin idxes in the first three return
    % values, i.e. < 0 or >= grid_info.nbin(d). In the last return value, these
    % invalid bin idxes are set to -1.
    %
    % <bin_idx> may be scalar (outputs are rows) or an array (outputs are [n_nbr, nb]).

    nbin = grid_info.nbin;
    offsets = grid_info.nbr_offsets;

    N_x_bins = nbin(1);
    N_y_bins = nbin(2);
    N_z_bins = nbin(3);

    b = double(bin_idx(:).');
    id_z = mod(b, N_z_bins);
    id_y = mod(floor(b / N_z_bins), N_y_bins);
    id_x = floor(b / (N_y_bins * N_z_bins));

    id_xs = offsets(1, :).' + id_x;   % [n_nbr, nb]
    id_ys = offsets(2, :).' + id_y;
    id_zs = offsets(3, :).' + id_z;

    binids = id_zs + id_ys * N_z_bins + id_xs * (N_y_bins * N_z_bins);
    binids(id_xs < 0 | id_xs >= N_x_bins | ...
           id_ys < 0 | id_ys >= N_y_bins | ...
           id_zs < 0 | id_zs >= N_z_bins) = -1;

    if isscalar(bin_idx)
        id_xs = id_xs.';
        id_ys = id_ys.';
        id_zs = id_zs.';
        binids = binids.';
    end
end
