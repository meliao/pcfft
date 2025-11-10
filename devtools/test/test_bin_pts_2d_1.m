addpath(genpath("../../pcfft"));

n_pts = 100000;
L = 2.0;
Lbd = [-1 -1 1 1];
% r points live on [-1, 1] x [-1, 1]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:);

% dx = 0.25, so the grid points are at
% x grid and y grid = [-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
dx = 0.25;
ngrid = [9 9];
% When we set nspread = 3, we expect
% x bins and y bins [-1.125, -0.375], [-0.375, 0.375], [0.375, 1.125]
nspread = 3;
nbin = [3 3];
[r_sorted, sorted_bin_ids, id_start] = bin_pts_2d(r, dx, ngrid, Lbd, nbin, nspread);
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
% for x = 2:size(id_start, 2)
%     idx = id_start(x);
%     before_idx = sorted_bin_ids(idx-1);
%     at_idx = sorted_bin_ids(idx);
%     assert(at_idx - before_idx == 1);
% end