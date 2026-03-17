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
    [idx_x, idx_y, binids] = intersecting_bins_2d(bin_idx, grid_info, proxy_info);

    % Once we remove invalid binids, there should be no repeats.
    valid_binids = binids(binids >= 0);
    assert(length(valid_binids) == length(unique(valid_binids)));

end

%% test_0b
% Larger test case with visualization. Same setup as test_SortInfo_2d_1


n_pts = 10000;
L = 2.0;
Lbd = [-1 1; 
        -1 1];
% r points live on [-1, 1] x [-1, 1]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;

% Get grid and proxy info
tol = 1e-6;
[grid_info, proxy_info] = get_grid(@log_kernel, ...
    struct('r', r), ...
    struct('r', r), ...
    tol, n_pts * 10);


% dx = 0.25, so the grid points are at
% dx = 0.25;
% ngrid = [9 9];
% % When we set nbinpts = 3, we expect
% % x bins and y bins [-1, -0.25], [-0.25, 0.5], [0.5, 1.]
% nbinpts = 3;
% nbin = [3 3];
% N_bin = nbin(1) * nbin(2);

% % Generate the GridInfo object. Need nbin, dx, Lbd, nspread, nbinpts, offset, dx 
% grid_info = GridInfo(Lbd, dx, 2*nbinpts + 1, nbinpts, dim, 0);
disp("test_intersecting_bins_2d: grid_info:");
disp(grid_info);
disp("test_intersecting_bins_2d: grid_info.nbin:");
disp(grid_info.nbin);

% Test the third return value is correct.

[~, ~, bin_4_intersecting_binids] = ...
    intersecting_bins_2d(4, grid_info, proxy_info);

% This should = [0 1 2 3 4 5 6 7 8
expected_bin_4_intersecting_binids = [0 1 2 3 4 5 6 7 8];
% disp("test_intersecting_bins_2d: For bin_idx 4, intersecting binids: ");
% disp(bin_4_intersecting_binids);
valid_bins = bin_4_intersecting_binids >= 0 & bin_4_intersecting_binids < N_bin;
valid_bins = bin_4_intersecting_binids(valid_bins);
% disp("test_intersecting_bins_2d: For bin_idx 4, valid intersecting binids: ");
% disp(valid_bins);
% assert(all(valid_bins == expected_bin_4_intersecting_binids));


[sort_info] = SortInfo(struct('r', r), grid_info.dx, grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);
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
N_bins = grid_info.nbin(1) * grid_info.nbin(2);
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
% In reality, we expect the unique sorted values of bin_0_intersecting_x 
% to be [-1, 0, 1]
expected_bin_0_intersecting = [-1 0 1];
unique_bin_0_intersecting_x = unique(bin_0_intersecting_x);
unique_bin_0_intersecting_y = unique(bin_0_intersecting_y);
assert(all(unique_bin_0_intersecting_x == expected_bin_0_intersecting));
assert(all(unique_bin_0_intersecting_y == expected_bin_0_intersecting));

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

assert(all(unique(bin_5_intersecting_x) == expected_bin_5_intersecting_x));
assert(all(unique(bin_5_intersecting_y) == expected_bin_5_intersecting_y));

close all;

%% test_0c
% Test that the returned bins for each query bin satisfy the proxy-circle
% intersection condition

rng(42);
n_src = 5000;
n_targ = 1500;
tol_0c = 1e-6;

src_info_0c = struct;
src_info_0c.r = rand(2, n_src) - 0.5;
src_info_0c.weights = rand(n_src, 1);

targ_info_0c = struct;
targ_info_0c.r = rand(2, n_targ) - 0.5;

[grid_info_0c, proxy_info_0c] = get_grid(@log_kernel, src_info_0c, targ_info_0c, tol_0c);

N_bin_0c = grid_info_0c.nbin(1) * grid_info_0c.nbin(2);
all_bins = 0:(N_bin_0c - 1);
tol_dist = 1e-10;

for bin_idx = 0:(N_bin_0c - 1)
    [~, ~, binids] = intersecting_bins_2d(bin_idx, grid_info_0c, proxy_info_0c);
    c1 = bin_center(bin_idx, grid_info_0c);
    valid_returned = binids(binids >= 0);
    assert(length(valid_returned) == length(unique(valid_returned)), ...
        sprintf('bin_idx %d: returned bins have repeats', bin_idx));

    % First returned bins must intersect
    for k = 1:length(valid_returned)
        c2 = bin_center(valid_returned(k), grid_info_0c);
        d = norm(c1 - c2);
        assert(d <= 2 * proxy_info_0c.radius + tol_dist, ...
            sprintf('bin %d -> bin %d dist=%.6f > 2*r=%.6f', ...
                    bin_idx, valid_returned(k), d, 2*proxy_info_0c.radius));
    end

    %  non-returned valid bins must not intersect
    non_returned = setdiff(all_bins, valid_returned);
    for k = 1:length(non_returned)
        c2 = bin_center(non_returned(k), grid_info_0c);
        d = norm(c1 - c2);
        assert(d > 2 * proxy_info_0c.radius - tol_dist, ...
            sprintf('bin %d -> bin %d dist=%.6f <= 2*r=%.6f but not returned', ...
                    bin_idx, non_returned(k), d, 2*proxy_info_0c.radius));
    end
end

