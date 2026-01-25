addpath(genpath("../../pcfft"));

n_pts = 100000;
L = 2.0;
Lbd = [-1 -1 1 1];
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

% Plot the x or y values of the sorted points to make sure the sorting is correct
% plot(r_sorted(2,:))


% Check the id_start array to make sure they indicate the correct indices.
% for x = 2:size(id_start, 2)
%     idx = id_start(x);
%     before_idx = sorted_bin_ids(idx-1);
%     at_idx = sorted_bin_ids(idx);
%     assert(at_idx - before_idx == 1);
% end