addpath(genpath("../../pcfft"));

% First check that the size is correct.
% bin_idx = 2;
% dx = 0.5;
% nbinpts = 3;
% nspread = 5;
% Lbd = [-2 2 
%        -2 2];
% Want spreading bins starting at [-2, -0.5, 1]
grid_info = struct;
grid_info.dx = 0.5;
grid_info.nbinpts = 3;
grid_info.nspread = 5;
grid_info.Lbd = [-2 2
                -2 2];
% From the notes and code in get_grid, we can derive
% the regular grid points.
bin_sidelen = grid_info.dx * grid_info.nbinpts;
n_bin = ceil(diff(grid_info.Lbd, 1, 2) / bin_sidelen);
disp("test: nbin: " + int2str(n_bin));
% Number of points padding each side
pad = ceil((grid_info.nspread - grid_info.nbinpts) / 2);
% Width below the bottom corner of Lbd to start the regular grid points
offset = pad * grid_info.dx - grid_info.dx / 2;
ngrid = n_bin * grid_info.nbinpts + pad * 2;
xx = grid_info.Lbd(1, 1) - offset + (0: ngrid(1) - 1) * grid_info.dx;
yy = grid_info.Lbd(2, 1) - offset + (0: ngrid(2) - 1) * grid_info.dx;


% Save some of these derived attributes in grid_info
grid_info.ngrid = ngrid;
grid_info.offset = offset;
grid_info.nbin = n_bin;

% Do a lot of points, then plot the gridpts on top of the 
% scattered points.
n_pts = 100000;
L = 4.0;
% r points live on [-2, 2]^2
rng(0);
r = (rand(2, n_pts) - 0.5) * L;

sort_info = SortInfo(r, grid_info.dx, grid_info.Lbd, ...
    grid_info.nbin, grid_info.nbinpts);
r_sorted = sort_info.r_srt;
sorted_bin_ids = sort_info.binid_srt;
id_start = sort_info.id_start;


% Get pts for a certain bin idx
bin_idx = 4;
[grid_pts, grid_ctr] = grid_pts_for_box_2d(bin_idx, grid_info);

disp("size of ctr")
disp(size(grid_ctr))

% Assert that grid_pts is the correct size.
assert(all(size(grid_pts) == [2, grid_info.nspread^2]));

[X, Y] = meshgrid(xx, yy);
% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
hold on;
scatter(grid_pts(1,:), grid_pts(2,:), 200, 'rx')
scatter(X(:), Y(:), 100, 'ko')
colormap('parula');
colorbar;