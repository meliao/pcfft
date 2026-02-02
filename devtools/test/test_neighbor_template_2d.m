% Makes sure get_addsub returns without error on a 2D input.
addpath(genpath('../../pcfft'));


rad = 10.0;
tol = 1e-10;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 13;
rng(0);
src_info_2d.r = (rand(2, n_src) - 0.5);
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
ntarg = 17;
targ_info_2d.r = rand(2, ntarg) - 0.5;


[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);

N_bin = grid_info.nbin(1) * grid_info.nbin(2);

% disp("test_intersecting_bins_2d: N_bin = " + int2str(N_bin));

s_temp = neighbor_template_2d(grid_info, proxy_info);



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
grid_info = GridInfo(Lbd, dx, nspread, nbinpts, dim);
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

[sort_info] = SortInfo(r, dx, Lbd, nbin, nbinpts);
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
temp_at_0 = neighbor_template_2d(grid_info, proxy_info);

% Compute the mean
disp("test_neighbor_template_2d: Mean of template at bin idx 0: ");
disp(mean(temp_at_0, 2));

% Now center at the center of bin idx 4.
ctr = bin_center(4, grid_info);
disp("test_neighbor_template_2d: Centering template at bin idx 4, ctr: ");
disp(ctr);
temp_at_4 = ctr + temp_at_0;

% Plot the regular grid points and plot the template points overtop them.
scatter(grid_info.r(1,:), grid_info.r(2,:), 10, 'b');

% Plot the template points
hold on;
scatter(temp_at_4(1,:), temp_at_4(2,:), 5, 'r', 'filled');
% close all;

% Confirm that each of the template points matches one of the grid points.
for i = 1:size(temp_at_4, 2)
    pt = temp_at_4(:, i);
    diffs = grid_info.r - pt;
    dists = sqrt(sum(diffs.^2, 1));
    min_dist = min(dists);
    assert(min_dist < 1e-12);
end
