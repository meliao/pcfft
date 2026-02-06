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
grid_info.rpad = pad;

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
bin_idx = 0;
[grid_pts, grid_ctr] = grid_pts_for_bin_2d(bin_idx, grid_info);

disp("size of ctr")
disp(size(grid_ctr))

% Assert that grid_pts is the correct size.
assert(all(size(grid_pts) == [2, grid_info.nbinpts^2]));

[X, Y] = meshgrid(xx, yy);
% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
hold on;
scatter(grid_pts(1,:), grid_pts(2,:), 200, 'rx')
scatter(X(:), Y(:), 100, 'ko')
colormap('parula');
colorbar;

%% test_1
% Check the row indices are computed correctly

% close all;

tol = 1e-7;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 1000;
rng(0);
src_info_2d.r = (rand(2, n_src) - 0.5) .^3;
src_info_2d.weights = rand(n_src, 1);


targ_info_2d = struct;
ntarg = 1700;
targ_info_2d.r = (rand(2, ntarg) - 0.5).^3;

[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);
disp("test: finished get_grid")
disp(grid_info.ngrid)

assert(grid_info.rpad > 1);

sort_info = SortInfo(src_info_2d.r, grid_info.dx, ...
    grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);
r_sorted = sort_info.r_srt;
sorted_bin_ids = sort_info.binid_srt;
id_start = sort_info.id_start;

close all;
disp("test: grid_info.Lbd");
disp(grid_info.Lbd);
scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
hold on;
scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
colorbar;
% Plot bin 1 with red x
[pts, center, row_indices] = grid_pts_for_bin_2d(1, grid_info);
scatter(pts(1,:), pts(2,:), 200, 'rx');
% Plot the associated spreading box
[boxpts, boxcenter, box_row_indices] = grid_pts_for_box_2d(1, grid_info);
scatter(boxpts(1,:), boxpts(2,:), 200, 'go');
close all;

N_x_bins = grid_info.nbin(1);
N_y_bins = grid_info.nbin(2);
N_bins = N_x_bins * N_y_bins;

disp("grid_info.ngrid")
disp(grid_info.ngrid)
% disp("id_start")
% disp(id_start)

% Assert that the row indices are computed correctly by checking 
% the points identified by src_info.r(:, row_indices) against pts.
for bin_id = 0:N_bins-1
    disp("test: On bin_id = " + int2str(bin_id))

    [pts, center, row_idxes] = grid_pts_for_bin_2d(bin_id, grid_info);
    % disp("test: Row_idxes: ")
    % disp(row_idxes)

    % pts_sliced should exactly match pts
    pts_sliced = grid_info.r(:, row_idxes);
    disp("test: pts_sliced size: " + int2str(size(pts_sliced)));
    disp(pts_sliced);
    disp("test: pts size: " + int2str(size(pts)));
    disp(pts);
    diffs = pts_sliced - pts;
    diffs = abs(diffs);
    disp("diffs")
    disp(diffs)
    
    % Scatter pts and pts_sliced
    % scatter(pts(1,:), pts(2,:), 100, 'ko');
    % hold on;
    % scatter(pts_sliced(1,:), pts_sliced(2,:), 200, 'rx');

    assert(all(size(pts_sliced) == size(pts)));
    assert(all(diffs(:) < 1e-15));
end