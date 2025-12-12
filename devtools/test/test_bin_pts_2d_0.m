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
% When we set nbinpts = 3, we expect
% x bins [-1.125, -0.375], [-0.375, 0.375], [0.375, 1.125]
% y bins [-0.625, 0.125], [0.125, 0.875]
% so we get nbin = [3 2]
nbin = [3 2];
N_bins = nbin(1) * nbin(2) + 1; % Total number of bins in this case.
nbinpts = 3;
[r_srt, binid_srt, ptid_srt, id_start] = bin_pts_2d(r, dx, Lbd, nbin, nbinpts);
c = 1:n_pts;
assert(all(size(r_srt) == size(r)));
assert(all(size(r, 2) == size(binid_srt, 2)));
disp(size(id_start))
assert(all(size(id_start, 2) == N_bins));


% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
% scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
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
[r_srt, binid_srt, ptid_srt, id_start] = bin_pts_2d(src_info.r, grid_info.dx, ...
 grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);


% N_x_bins = ceil(grid_info.ngrid(1)/nbin);
% N_y_bins = ceil(grid_info.ngrid(2)/nbin);
% N_bins = N_x_bins * N_y_bins;

% assert(all(size(id_start) == [1 N_bins]));