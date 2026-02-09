% Plot a failing test cases I've seen in the wild
addpath(genpath('../../pcfft'));


% tol = 1e-10;
% dim = 2;

% k = @(s,t) log_kernel(s,t);

% src_info_2d = struct;
% n_src = 2;
% rng(0);
% src_info_2d.r = (rand(2, n_src) - 0.5);
% src_info_2d.weights = zeros(n_src, 1);
% src_info_2d.weights(2) = 1.0;

% targ_info_2d = struct;
% ntarg = 17;
% targ_info_2d.r = src_info_2d.r;


% [grid_info, proxy_info] = get_grid(@log_kernel, ...
%     src_info_2d, targ_info_2d, tol);

% sort_info = SortInfo(src_info_2d.r, grid_info.dx, ...
%     grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);
% r_sorted = sort_info.r_srt;
% sorted_bin_ids = sort_info.binid_srt;
% id_start = sort_info.id_start;

% disp("sorted_bin_ids")
% disp(sorted_bin_ids)
% scatter(grid_info.r(1, :), grid_info.r(2,:), 'k.');
% hold on;
% scatter(r_sorted(1,:), r_sorted(2,:), 20, sorted_bin_ids, 'filled');
% colorbar;

% for bin_id = sorted_bin_ids
%     [pts_i, center, row_idxes] = grid_pts_for_box_2d(bin_id, grid_info, 1);
%     scatter(pts_i(1, :), pts_i(2, :), 'ro');
% end

%% Second case: 2D problem with 2 points.
rng(3);
n_src = 2;
half_side_len = 1.0;
source_pts = (rand(2,n_src) - 0.5) * half_side_len;
% Cube the uniform IID points so there is a non-uniform distribution
source_pts(1,:) = source_pts(1,:).^3;
source_pts(2,:) = source_pts(2,:).^3;
disp("source_pts:")
disp(source_pts)

src_weights = rand(n_src,1);
src_weights = src_weights(:);

src_info = struct;
src_info.r = source_pts;

% Make the "target points" the same as the source points.
targ_info = struct;
targ_info.r = source_pts;

n_nbr = 500;
k = @(s,t) log_kernel(s,t);
tol = 1e-05;

[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);

sort_info = SortInfo(src_info.r, grid_info.dx, ...
    grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);
r_sorted = sort_info.r_srt;
sorted_bin_ids = sort_info.binid_srt;
id_start = sort_info.id_start;
disp("test: sorted_bin_ids")
disp(sorted_bin_ids)

disp("test: id_start")
disp(id_start)

for bin_id = sorted_bin_ids
    [pts_i, center, row_idxes] = grid_pts_for_box_2d(bin_id, grid_info);
    % pts_sliced should exactly match pts
    pts_sliced = grid_info.r(:, row_idxes);
    disp("test: pts_sliced");
    disp(pts_sliced);
    disp("test: pts_i");
    disp(pts_i);
    diffs = pts_sliced - pts_i;
    % disp("diffs")
    % disp(diffs)
    assert(all(size(pts_sliced) == size(pts_i)));
    assert(all(diffs(:) < 1e-15));

    % Assert that all of the row_idxes are < = n_grid_pts
    n_grid_pts = size(grid_info.r, 2);
    assert(all(row_idxes <= n_grid_pts));
end