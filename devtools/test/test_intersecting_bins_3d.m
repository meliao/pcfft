% Test that intersecting_bins_3d returns a valid neighborhood for every bin.
addpath(genpath('../../pcfft'));

dx = 0.25;
nbinpts = 3;
dim = 3;

Lbd = [-1  1
       -1  1
       -1  1];

% Spoof a proxy_info with a reasonable radius.
proxy_info = struct;
proxy_info.radius = 0.5;

grid_info = GridInfo(Lbd, dx, 2*nbinpts + 1, nbinpts, dim, 0, proxy_info.radius);



N_bin = grid_info.nbin(1) * grid_info.nbin(2) * grid_info.nbin(3);

for bin_idx = 0:(N_bin - 1)
    [~, ~, ~, binids] = intersecting_bins_3d(bin_idx, grid_info);

    valid_binids = binids(binids >= 0);
    assert(length(valid_binids) <= N_bin, ...
        sprintf('bin_idx %d: neighborhood size %d exceeds N_bin %d', ...
                bin_idx, length(valid_binids), N_bin));
end

