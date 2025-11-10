addpath(genpath("../../pcfft"));

% First check that the size is correct.
% bin_idx = 2;
% dx = 0.5;
% ngrid = [24 24];
% Lbd = [0 0 1 1];
% nbin = 7;
% pts = grid_pts_for_bin_2d(bin_idx, dx, ngrid, Lbd, nbin);
% assert(all(size(pts) == [2 nbin^2]));

% Do a lot of points, then plot the gridpts on top of the 
% scattered points.
n_pts = 100000;
L = 2.0;
% Lbd = [-1 -0.5 1 0.5];
Lbd = [-1 1
 -1 1];
% r points live on [-1, 1]^2
rng(0);
r = (rand(2, n_pts) - 0.5) * L;

% dx = 0.25, so the grid points are at
% grid = [-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
dx = 0.25;
% ngrid = [9 5];
ngrid = [9 9];
% When we set nspread = 3, we expect
% bins [-1.125, -0.375], [-0.375, 0.375], [0.375, 1.125]
nspread = 3;
nbin = [3 3];
[r_sorted, sorted_bin_ids, id_start] = bin_pts_2d(r, dx, ngrid, Lbd, nbin, nspread);

grid_info = struct;
grid_info.dx = dx;
grid_info.ngrid = ngrid;
grid_info.rpad = 2;
grid_info.nbin = nbin;
grid_info.nspread = nspread;
grid_info.Lbd = Lbd;


% Get pts for a certain bin idx
bin_idx = 4;
[grid_pts, grid_ctr] = grid_pts_for_bin_2d(bin_idx, grid_info);

disp("size of ctr")
disp(size(grid_ctr))

xgrid = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
[X, Y] = meshgrid(xgrid, xgrid);
% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
hold on;
scatter(grid_pts(1,:), grid_pts(2,:), 100, 'rx')
scatter(X(:), Y(:), 100, 'ko')
colormap('parula');
colorbar;