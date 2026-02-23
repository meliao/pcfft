addpath(genpath("../../pcfft"));


% Do a lot of points, then plot the gridpts on top of the 
% scattered points.
n_pts = 100000;
L = 2.0;
Lbd = [-1 1
    -0.5 0.5];
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:)/2;

% dx = 0.25, so the spreading bins are
% x bins [-1, -0.25], [-0.25, 0.5], [0.5, 1]
% y bins [-0.5, 0.25], [0.25, 0.5]
dx = 0.25;
nbin = [3 2];
nspread = 6;
nbinpts = 3;

grid_info = struct;
grid_info.dx = dx;
grid_info.nbinpts = nbinpts;
grid_info.nspread = nspread;
bin_sidelen = dx * grid_info.nbinpts;

% Number of spreading bins in each dimension.
n_bin = ceil(diff(Lbd, 1, 2) / bin_sidelen);
% Number of points padding each side
pad = ceil((grid_info.nspread - grid_info.nbinpts) / 2);
% Width below the bottom corner of Lbd to start the regular grid points
offset = pad * dx - dx / 2;

ngrid = n_bin * grid_info.nbinpts + pad * 2;

grid_info.ngrid = ngrid;
grid_info.rpad = 2;
grid_info.Lbd = Lbd;
grid_info.nbin = nbin;
grid_info.nspread = nspread;
grid_info.nbinpts = nbinpts;
pad = ceil((grid_info.nspread - grid_info.nbinpts) / 2);
% Width below the bottom corner of Lbd to start the regular grid points
offset = pad * grid_info.dx - grid_info.dx / 2;
grid_info.offset = offset;

sort_info = SortInfo(struct('r', r), grid_info.dx, grid_info.Lbd, nbin, nbinpts);
r_sorted = sort_info.r_srt;
sorted_bin_ids = sort_info.binid_srt;
id_start = sort_info.id_start;




% Get pts for a certain bin idx
bin_idx = 4;
[grid_pts, grid_ctr] = grid_pts_for_box_2d(bin_idx, grid_info);

disp("size of ctr")
disp(size(grid_ctr))

xgrid = Lbd(1, 1) - offset + (0: ngrid(1) - 1) * dx;
ygrid = Lbd(2, 1) - offset + (0: ngrid(2) - 1) * dx;
[X, Y] = meshgrid(xgrid, ygrid);
% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
hold on;
scatter(grid_pts(1,:), grid_pts(2,:), 100, 'rx')
scatter(X(:), Y(:), 100, 'ko')
colormap('parula');
colorbar;