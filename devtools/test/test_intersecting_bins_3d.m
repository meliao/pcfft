% Test that intersecting_bins_3d returns a valid neighborhood for every bin.
addpath(genpath('../../pcfft'));

dx = 0.25;
nbinpts = 3;
dim = 3;

Lbd = [-1  1
       -1  1
       -1  1];

grid_info = GridInfo(Lbd, dx, 2*nbinpts + 1, nbinpts, dim, 0);

% Spoof a proxy_info with a reasonable radius.
proxy_info = struct;
proxy_info.radius = 0.5;

N_bin = grid_info.nbin(1) * grid_info.nbin(2) * grid_info.nbin(3);

for bin_idx = 0:(N_bin - 1)
    [~, ~, ~, binids] = intersecting_bins_3d(bin_idx, grid_info, proxy_info);
    assert(length(binids) <= N_bin, ...
        sprintf('bin_idx %d: neighborhood size %d exceeds N_bin %d', ...
                bin_idx, length(binids), N_bin));
end

disp('test_intersecting_bins_3d: PASSED');
