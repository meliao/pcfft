addpath(genpath("../../pcfft"));

% First check that the size is correct.
bin_idx = 2;
dx = 0.5;
ngrid = [24 24];
Lbd = [0 0 1 1];
nbin = 7;
pts = grid_pts_for_bin_2d(bin_idx, dx, ngrid, Lbd, nbin);
assert(all(size(pts) == [2 nbin^2]));

% Do a lot of points, then plot the gridpts on top of the 
% scattered points.
n_pts = 10000;
L = 2.0;
Lbd = [-1 -0.5 1 0.5];
% r points live on [-1, 1] x [-0.5 0.5]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:) / 2;
nbin = 2;

dx = 0.5;
ngrid = [5 3];
[r_sorted, sorted_bin_ids, id_start] = bin_pts_2d(r, dx, ngrid, Lbd);

% Get pts for a certain bin idx
bin_idx = 0;
grid_pts = grid_pts_for_bin_2d(bin_idx, dx, ngrid, Lbd, nbin);

% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
hold on;
scatter(grid_pts(1,:), grid_pts(2,:), 100, 'rx')
colormap('parula');
colorbar;