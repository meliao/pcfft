addpath(genpath("../../pcfft"));

n_pts = 100000;
L = 2.0;
Lbd = [-1 -0.5 1 0.5];
% r points live on [-1, 1] x [-0.5 0.5]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:) / 2;

% dx = 0.25, so the grid points are at
% x grid = [-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
% y grid = [-0.5, -0.25, 0.0, 0.25 0.5]
dx = 0.25;
ngrid = [9 5];
% When we set nbin = 3, we expect
% x bins [-1.125, -0.375], [-0.375, 0.375], [0.375, 1.125]
% y bins [-0.625, 0.125], [0.125, 0.875]
nbin = 3;
[r_sorted, sorted_bin_ids, id_start] = bin_pts_2d(r, dx, ngrid, Lbd, nbin);
c = 1:n_pts;
assert(all(size(r_sorted) == size(r)));
assert(all(size(r, 2) == size(sorted_bin_ids, 2)));


% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
colormap('parula');
colorbar;

% Plot the x or y values of the sorted points to make sure the sorting is correct
% plot(r_sorted(2,:))


% Check the id_start array to make sure they indicate the correct indices.
for x = 2:size(id_start, 2)
    idx = id_start(x);
    before_idx = sorted_bin_ids(idx-1);
    at_idx = sorted_bin_ids(idx);
    assert(at_idx - before_idx == 1);
end