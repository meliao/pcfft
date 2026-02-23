% Makes sure get_addsub returns without error on a 2D input.
addpath(genpath('../../pcfft'));

close all;

tol = 1e-10;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 137;
rng(0);
src_info_2d.r = (rand(2, n_src) - 0.5);
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
ntarg = 173;
targ_info_2d.r = rand(2, ntarg) - 0.5;


[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);

N_bin = grid_info.nbin(1) * grid_info.nbin(2);

% disp("test_intersecting_bins_2d: N_bin = " + int2str(N_bin));

[nbr_binids, nbr_gridpts, nbr_grididxes, bin_idx] = neighbor_template_2d(grid_info, proxy_info, 8);

assert( bin_idx == 8);

n_x = sqrt(length(nbr_binids));
expected_n_pts = n_x * grid_info.nbinpts + 2 * grid_info.rpad;
disp("test_neighbor_template_2d: expected_n_pts: ");
disp(expected_n_pts^2);
disp("test_neighbor_template_2d: nbr_gridpts size: ");
disp(size(nbr_gridpts));
disp("test_neighbor_template_2d: nbr_grididxes size: ");
disp(size(nbr_grididxes));
assert(size(nbr_gridpts, 2) == expected_n_pts^2);
assert(size(nbr_grididxes, 2) == expected_n_pts^2);



%% test_0b
% Larger test case with visualization. Same setup as test_SortInfo_2d_1


n_pts = 100000;
L = 2.0;
Lbd = [-1 1; -1 1];
% r points live on [-1, 1] x [-1, 1]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:);

% src_pts = struct('r', r);

% [grid_info, proxy_info] = get_grid(@log_kernel, ...
%     src_pts, src_pts, tol);

% dx = 0.25, so the grid points are at
dx = 0.25;
ngrid = [9 9];
% When we set nbinpts = 3, we expect
% x bins and y bins [-1, -0.25], [-0.25, 0.5], [0.5, 1.]
nbinpts = 3;
nbin = [3 3];
nspread = 2 * nbinpts ;

% Create a GridInfo object. Need nbin, dx, Lbd, nspread, nbinpts, offset, dx, ngrid
grid_info = GridInfo(Lbd, dx, nspread, nbinpts, dim, -1);
% grid_info = struct;
% grid_info.nbin = nbin;
% grid_info.dx = dx;
% grid_info.Lbd = Lbd;
% grid_info.nspread = 2*nbinpts + 1;
% grid_info.nbinpts = nbinpts;
% pad = ceil((grid_info.nspread - nbinpts)/2);
% grid_info.offset = pad * dx - dx/2;
% grid_info.dx = dx;

% spoof the ProxyInfo object. Need radius
proxy_info = struct;
proxy_info.radius = 0.5;

[sort_info] = SortInfo(struct('r', r), dx, Lbd, nbin, nbinpts);
r_srt = sort_info.r_srt;
binid_srt = sort_info.binid_srt;
ptid_srt = sort_info.ptid_srt;
id_start = sort_info.id_start;
c = 1:n_pts;
assert(all(size(r_srt) == size(r)));
assert(all(size(r, 2) == size(binid_srt, 2)));


% Plot the sorted points and color by the bin
% to make sure the bin assignment looks correct
% scatter(r_srt(1,:), r_srt(2,:), 20, binid_srt, 'filled');
% colormap('parula');
% colorbar;

% Draw an x at the center of each bin
% hold on;
% N_bins = nbin(1) * nbin(2);
% for bin_idx = 0:(N_bins - 1)
%     center = bin_center(bin_idx, grid_info);
%     scatter(center(1), center(2), 100, 'x');

%     % Draw the proxy ring
%     proxypts = get_ring_points(100, proxy_info.radius, center);
%     plot(proxypts(1,:), proxypts(2,:), 'k-');
% end

% Figure shows that bin idx 0 only intersects with 0, 1, 3. 
[nbr_binids, nbr_gridpts, nbr_grididxes, bin_idx] = neighbor_template_2d(grid_info, proxy_info, 0);

% Now center at the center of bin idx 0.
ctr = bin_center(bin_idx, grid_info);
disp("test_neighbor_template_2d: Centering template at bin idx 4, ctr: ");
disp(ctr);
temp_at_4 = nbr_gridpts;

% Plot the regular grid points and plot the template points overtop them.
scatter(grid_info.r(1,:), grid_info.r(2,:), 10, 'b');

% Plot the template points
hold on;
scatter(temp_at_4(1,:), temp_at_4(2,:), 5, 'r', 'filled');


% Plot text lables for the valid nbr idxes
text(temp_at_4(1,:), temp_at_4(2,:), string(nbr_grididxes), 'Color', 'k');

% Confirm that each of the template points matches one of the grid points, after
% filtering
oob_idxes = nbr_grididxes > grid_info.ngrid(1) * grid_info.ngrid(2);
valid_temp_pts = temp_at_4(:, ~oob_idxes);

% Plot the valid template points in green
scatter(valid_temp_pts(1,:), valid_temp_pts(2,:), 20, 'g', 'filled');
% close all;

for i = 1:size(valid_temp_pts, 2)
    pt = valid_temp_pts(:, i);
    diffs = grid_info.r - pt;
    dists = sqrt(sum(diffs.^2, 1));
    min_dist = min(dists);
    assert(min_dist < 1e-12);
end


% Check to make sure that the nbr_grididxes for valid points are correct.
valid_nbr_idxes = nbr_grididxes(~oob_idxes);

% Plot text lables for the valid nbr idxes
text(valid_temp_pts(1,:), valid_temp_pts(2,:), string(valid_nbr_idxes), 'Color', 'k');

% Plot the spreading box for bin idx 0
[boxpts, boxctr, boxidxes] = grid_pts_for_box_2d(0, grid_info);
disp("test_neighbor_template_2d: boxidxes: ");
disp(boxidxes);
boxpts2 = grid_info.r(:, boxidxes);
% Check equality between boxpts and boxpts2
diffs = boxpts - boxpts2;
dists = sqrt(sum(diffs.^2, 1));
assert(max(dists) < 1e-12);
scatter(boxpts2(1,:), boxpts2(2,:), 20, 'm', 'filled');


for i = 1:size(valid_temp_pts, 2)
    pt = valid_temp_pts(:, i);
    idx = valid_nbr_idxes(i);
    % disp("test_neighbor_template_2d: Checking pt " + int2str(i) + " at idx " + int2str(idx));
    % disp("pt: ");
    % disp(pt);
    % disp("grid pt: ");
    grid_pt = grid_info.r(:, idx);
    % disp(grid_pt);
    dist = norm(pt - grid_pt);
    assert(dist < 1e-12);
end


% close all;
%% test_0c

% Check that for a given A_spread_s matrix, indexing it with the 
% nbr_grididxes gets all of the relevant entries for a given bin_idx.
rad = 2.0;

% Set up two source points
rng(4);
n_src = 2;
source_pts = [-0.0531 0.00223
                -0.00497 0.4143];

% target points are close but not exactly = source points
target_pts = [-0.0491 0.00234 0.1
                -0.00678 0.4102 0.1];
disp(size(target_pts));


src_weights = zeros(n_src,1);
src_weights(1) = 1.0;
src_weights = src_weights(:);

% Define the kernel
k = @(s,t) log_kernel(s,t);

K_src_to_target = log_kernel(struct('r',source_pts), struct('r',target_pts));

target_vals = K_src_to_target * src_weights;
n_nbr = 3; % 10000 points / 500 is approximately 20 boxes

src_info = struct;
src_info.r = source_pts;
targ_info = struct;
targ_info.r = target_pts;
tol = 1e-08;

[grid_info, proxy_info] = get_grid(k, src_info, targ_info, tol, n_nbr);

[A_spread_s, K_src_to_reg, sort_info_s] = get_spread(k, k, src_info, ...
    grid_info, proxy_info);

[A_spread_t, K_targ_to_reg, sort_info_t] = get_spread(k, k, targ_info, ...
    grid_info, proxy_info);

n_bins = grid_info.nbin(1) * grid_info.nbin(2);

for bin_idx = 0:(n_bins - 1)
    [nbr_binids, nbr_gridpts, nbr_grididxes, ~] = neighbor_template_2d(grid_info, proxy_info, bin_idx);
    disp("test_neighbor_template_2d: Checking bin idx " + int2str(bin_idx));
    disp("nbr_grididxes: ");
    disp(nbr_grididxes);

    % Filter out the nbr_gridpts which correspond to invalid nbr_grididxes
    oob_idxes = nbr_grididxes > grid_info.ngrid(1) * grid_info.ngrid(2);
    valid_nbr_grididxes = nbr_grididxes(~oob_idxes);
    valid_nbr_gridpts = nbr_gridpts(:, ~oob_idxes);
    
    % Make a new figure with the grid points and the template points for this bin idx
    figure;
    scatter(grid_info.r(1,:), grid_info.r(2,:), 10, 'b');
    hold on;
    scatter(valid_nbr_gridpts(1,:), valid_nbr_gridpts(2,:), 20, 'rx');
    title("Bin idx " + int2str(bin_idx));
    xlabel("x");
    ylabel("y");
end

close all;