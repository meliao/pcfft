addpath(genpath("../../pcfft"));

n_pts = 100000;
L = 2.0;
Lbd = [-1 -0.5 1 0.5];
% r points live on [-1, 1] x [-0.5 0.5]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:) / 2;

% dx = 0.25, so the grid points are at
dx = 0.25;
ngrid = [9 5];
% When we set nbinpts = 3, we expect
% x bins [-1, -0.25], [-0.25, 0.5], [0.5, 1.0]
% y bins [-0.5, 0.25], [0.25, 0.5]
% so we get nbin = [3 2]
nbin = [3 2];
N_bins = nbin(1) * nbin(2) + 1; % Total number of bins in this case.
nbinpts = 3;
sort_info = SortInfo(struct('r', r), dx, Lbd, nbin, nbinpts);
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
% scatter(r_srt(1,:), r_srt(2,:), 20, binid_srt, 'filled');
% colormap('parula');
% colorbar;

% Plot the x or y values of the sorted points to make sure the sorting is correct
% plot(r_sorted(2,:))


% Check the id_start array to make sure they 
% indicate the correct indices. This grid is super dense so there
% will not be any empty bins.
% disp(id_start);
% for x = 2:size(id_start, 2)
%     idx = id_start(x);
%     before_idx = sorted_bin_ids(idx-1);
%     at_idx = sorted_bin_ids(idx);
%     assert(at_idx - before_idx == 1);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

rng(0);
src_info = struct;
n_src = 13;
src_info.r = rand(2, n_src);
tol = 1e-10;

targ_info = struct;
n_targ = 17;
targ_info.r = rand(2, n_targ);

% Get a realistic grid which has empty bins
[grid_info, proxy_info] = get_grid(@log_kernel, src_info, targ_info, tol);
sort_info = SortInfo(src_info, grid_info.dx, ...
 grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);

r_srt = sort_info.r_srt;
binid_srt = sort_info.binid_srt;
ptid_srt = sort_info.ptid_srt;
id_start = sort_info.id_start;
% N_x_bins = ceil(grid_info.ngrid(1)/nbin);
% N_y_bins = ceil(grid_info.ngrid(2)/nbin);
% N_bins = N_x_bins * N_y_bins;

% assert(all(size(id_start) == [1 N_bins]));