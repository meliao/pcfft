% Makes sure get_addsub returns without error on a 2D input.
addpath(genpath('../../pcfft'));


rad = 10.0;
tol = 1e-10;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 13;
rng(0);
src_info_2d.r = (rand(2, n_src) - 0.5);
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
ntarg = 17;
targ_info_2d.r = rand(2, ntarg) - 0.5;


[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);

N_bin = grid_info.nbin(1) * grid_info.nbin(2);

% disp("test_intersecting_bins_2d: N_bin = " + int2str(N_bin));

% valid bin_idxes should be between 0 and N_bin - 1
for bin_idx = 0:(N_bin - 1)
    [idx_x, idx_y] = intersecting_bins_2d(bin_idx, grid_info, proxy_info);
    % disp("test_intersecting_bins_2d: For bin_idx " + int2str(bin_idx) + ...
    %     ", intersecting bins: ");
    % disp(bin_idxes);


end

%% test_0b
% Larger test case with visualization. Same setup as test_SortInfo_2d_1


n_pts = 100000;
L = 2.0;
Lbd = [-1 1; 
        -1 1];
% r points live on [-1, 1] x [-1, 1]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:);

% dx = 0.25, so the grid points are at
dx = 0.25;
ngrid = [9 9];
% When we set nbinpts = 3, we expect
% x bins and y bins [-1, -0.25], [-0.25, 0.5], [0.5, 1.]
nbinpts = 3;
nbin = [3 3];

% Generate the GridInfo object. Need nbin, dx, Lbd, nspread, nbinpts, offset, dx 
grid_info = GridInfo(Lbd, dx, 2*nbinpts + 1, nbinpts, dim);
disp("test_intersecting_bins_2d: grid_info:");
disp(grid_info);

% Test the third return value is correct.

[~, ~, bin_4_intersecting_binids] = ...
    intersecting_bins_2d(4, grid_info, proxy_info);

% This should = [0 1 2 3 4 5 6 7 8
expected_bin_4_intersecting_binids = [0 1 2 3 4 5 6 7 8];
disp("test_intersecting_bins_2d: For bin_idx 4, intersecting binids: ");
disp(bin_4_intersecting_binids);
assert(all(bin_4_intersecting_binids == expected_bin_4_intersecting_binids));

% spoof the ProxyInfo object. Need radius
proxy_info = struct;
proxy_info.radius = 0.5;

[sort_info] = SortInfo(r, dx, Lbd, nbin, nbinpts);
r_srt = sort_info.r_srt;
binid_srt = sort_info.binid_srt;
ptid_srt = sort_info.ptid_srt;
id_start = sort_info.id_start;
c = 1:n_pts;
assert(all(size(r_srt) == size(r)));
assert(all(size(r, 2) == size(binid_srt, 2)));


% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter(r_srt(1,:), r_srt(2,:), 20, binid_srt, 'filled');
colormap('parula');
colorbar;

% Draw an x at the center of each bin
hold on;
N_bins = nbin(1) * nbin(2);
for bin_idx = 0:(N_bins - 1)
    center = bin_center(bin_idx, grid_info);
    scatter(center(1), center(2), 100, 'x');

    % Draw the proxy ring
    proxypts = get_ring_points(100, proxy_info.radius, center);
    plot(proxypts(1,:), proxypts(2,:), 'k-');
end

% Figure shows that bin idx 0 only intersects with 0, 1, 3.
[bin_0_intersecting_x, bin_0_intersecting_y] = intersecting_bins_2d(0, grid_info, proxy_info);

% Expect bin_0_intersecting_x = [-1, 0, 1]
disp("test_intersecting_bins_2d: For bin_idx 0, intersecting bins x: ");
disp(bin_0_intersecting_x);
% expected_bin_0_intersecting_y = [-1, 0, 1];
disp("test_intersecting_bins_2d: For bin_idx 0, intersecting bins y: ");
disp(bin_0_intersecting_y);
expected_bin_0_intersecting = [-1 0 1];
assert(all(bin_0_intersecting_x == expected_bin_0_intersecting));
assert(all(bin_0_intersecting_y == expected_bin_0_intersecting));

% Figure shows that bin idx 5 intersects with 2, 4, 5, 8.
% bin_idx 5 corresponds to (id_x, id_y) = (1, 2)
% So we expect intersecting bins to be
% id_x in [0, 1, 2]
% id_y in [1, 2, 3]
[bin_5_intersecting_x, bin_5_intersecting_y] = intersecting_bins_2d(5, grid_info, proxy_info);
disp("test_intersecting_bins_2d: For bin_idx 4, intersecting bins: ");
disp(bin_5_intersecting_x);
disp(bin_5_intersecting_y);
expected_bin_5_intersecting_x = [0 1 2];
expected_bin_5_intersecting_y = [1 2 3];

assert(all(bin_5_intersecting_x == expected_bin_5_intersecting_x));
assert(all(bin_5_intersecting_y == expected_bin_5_intersecting_y));

close all;

%% test_0c

