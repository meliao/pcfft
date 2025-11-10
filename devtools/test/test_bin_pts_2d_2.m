% Saw this error when running an accuracy test.
% This tests that id_start has positive entries
addpath(genpath("../../pcfft"));

% Set up a random source
rng(3);
n_src = 10000;
half_side_len = 1.0;
source_pts = (rand(2,n_src) - 0.5) * half_side_len;
% Cube the uniform IID points so there is a non-uniform distribution
source_pts(1,:) = source_pts(1,:).^3;
source_pts(2,:) = source_pts(2,:).^3;

src_weights = rand(n_src,1);
src_weights = src_weights(:);

src_info = struct;
src_info.r = source_pts;

% Make the "target points" the same as the source points.
targ_info = struct;
targ_info.r = source_pts;

n_nbr = 500; % 10000 points / 500 is approximately 20 boxes
tol = 1e-04;
k = @(s,t) log_kernel(s,t);
[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);
nbin = grid_info.nbin;

[r_sorted, bin_idxes, id_start] = bin_pts_2d(src_info.r, ...
    grid_info.dx, grid_info.ngrid, grid_info.Lbd, grid_info.nbin, grid_info.nspread);

% Display the last 10 entries of id_start
disp("id_start:")
disp(id_start(end-10:end))

% Display the last 10 entries of bin_idxes
disp("bin_idxes:")
disp(bin_idxes(end-10:end))

% Assert that all of the id_start are positive
assert(all(id_start > 0));

% Assert that the last id_start is equal to n_src + 1
assert(id_start(end) == n_src + 1);

% Assert that the bin_idxes are all within the correct
% range.

max_bin_idx = grid_info.nbin(2) * grid_info.nbin(1) - 1;
assert(all(bin_idxes <= max_bin_idx));