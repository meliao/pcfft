addpath(genpath("../../pcfft"));

n_pts = 10000;
L = 2.0;
Lbd = [-1 1 % xmin xmax
       -0.5 0.5 % ymin ymax
       1 2.0]; % zmin zmax
% r points live on [-1, 1] x [-0.5 0.5] x [1, 2.0]
rng(0);
r = (rand(3, n_pts) - 0.5) * L;
r(2,:) = r(2,:) / 2;
r(3,:) = r(3,:) / 2 + 1;

% dx = 0.25, so the grid points are at
dx = 0.25;
ngrid = [9 5 5];
% When we set nbinpts = 3, we expect
% x bins [-1, -0.25], [-0.25, 0.5], [0.5, 1.0]
% y bins [-0.5, 0.25], [0.25, 0.5]
% z bins [1.0, 1.75], [1.75, 2.0]
% so we get nbin = [3 2 2]
nbin = [3 2 2];
N_bins = nbin(1) * nbin(2) * nbin(3) + 1;
nbinpts = 3;
src_info = struct('r', r);
sort_info = SortInfo(src_info, dx, Lbd, nbin, nbinpts);
r_srt = sort_info.r_srt;
binid_srt = sort_info.binid_srt;
id_start = sort_info.id_start;


c = 1:n_pts;
assert(all(size(r_srt) == size(r)));
assert(all(size(r, 2) == size(binid_srt, 2)));
disp(size(id_start))
assert(all(size(id_start, 2) == N_bins));

% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter3(r_srt(1,:), r_srt(2,:), r_srt(3,:), 20, binid_srt, 'filled');
colormap('parula');
colorbar;
