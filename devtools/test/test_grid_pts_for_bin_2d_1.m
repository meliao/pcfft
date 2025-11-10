addpath(genpath("../../pcfft"));

% Set up a geometry with ~ 100 bins

tol = 1e-3;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 1000;
rng(0);
src_info_2d.r = (rand(2, n_src) - 0.5) .^3;
src_info_2d.weights = rand(n_src, 1);


targ_info_2d = struct;
ntarg = 17;
targ_info_2d.r = (rand(2, ntarg) - 0.5).^3;

[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);
disp("test: finished get_grid")
disp(grid_info.ngrid)

[r_sorted, sorted_bin_ids, id_start] = bin_pts_2d(src_info_2d.r, ...
 grid_info.dx, grid_info.ngrid, grid_info.Lbd, grid_info.nbin, grid_info.nspread);

disp("test: grid_info.Lbd");
disp(grid_info.Lbd);
% scatter(grid_info.r(1,:), grid_info.r(2,:), 'k.');
% hold on;
% scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
% colorbar;

disp("test: finished bin_pts_2d")

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
    % disp("test: pts_sliced");
    % disp(pts_sliced);
    % disp("test: pts");
    % disp(pts);
    diffs = pts_sliced - pts;
    % disp("diffs")
    % disp(diffs)
    assert(all(size(pts_sliced) == size(pts)));
    assert(all(diffs(:) < 1e-15));
end