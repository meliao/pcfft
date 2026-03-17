% Tests for abstract_neighbor_spreading_2D.
addpath(genpath('../../pcfft'));

%% test_0: 
% Basic size and centering checks on a hand-constructed grid.
%
% dx=0.5, nbinpts=2, nspread=4, so each spreading box is 4x4=16 pts.

dim = 2;
Lbd = [-1 1; -1 1];
dx = 0.5;
nbinpts = 2;
nspread = 4;
grid_info = GridInfo(Lbd, dx, nspread, nbinpts, dim, -1);

proxy_info = struct;
proxy_info.radius = 0.5;

[box_pts, spreading_template_pts] = abstract_neighbor_spreading_2D(grid_info, proxy_info);

% box_pts must be [2, nspread^2]
assert(all(size(box_pts) == [2, nspread^2]), ...
    'box_pts must be [2, nspread^2]');

% box_pts must be centered at center_bin's spreading box center
[~, expected_box_center] = grid_pts_for_box_2d(grid_info.center_bin, grid_info);
box_center = mean(box_pts, 2);
assert(norm(box_center - expected_box_center) < 1e-10, ...
    'box_pts must be centered at center_bin spreading box center');

% No duplicate points in spreading_template_pts
[~, uid] = unique(spreading_template_pts.', 'rows');
assert(length(uid) == size(spreading_template_pts, 2), ...
    'spreading_template_pts must have no duplicate points');

% box_pts must be a subset of spreading_template_pts (central bin is included)
% spreading_template_pts is relative to box_center, so compare box_pts - box_center
for i = 1:size(box_pts, 2)
    dists = sqrt(sum((spreading_template_pts - (box_pts(:, i) - box_center)).^2, 1));
    assert(min(dists) < 1e-12, ...
        sprintf('box_pts column %d not found in spreading_template_pts', i));
end

% spreading_template_pts must have at least nspread^2 points
assert(size(spreading_template_pts, 2) >= nspread^2, ...
    'spreading_template_pts must have at least nspread^2 points');

% spreading_template_pts must be separated by about dx
dists = sqrt(sum(diff(spreading_template_pts, 1, 2).^2, 1));
min_dist = min(dists);
assert(min_dist + 1e-12 >= dx, ...
    sprintf('Minimum distance between spreading_template_pts is %g, less than dx=%g', min_dist, dx));


%% test_1: 
% Template size matches an interior bin from neighbor_template_2d.
%
% For a grid large enough that interior bins exist, the number of valid
% (in-bounds) points returned by neighbor_template_2d for an interior bin
% must equal size(spreading_template_pts, 2).

rng(0);
n_src = 2000;
tol = 1e-8;

scale = 4.0;

src_info = struct;
src_info.r = (rand(2, n_src) - 0.5) * scale;


[grid_info_1, proxy_info_1] = get_grid(@log_kernel, src_info, src_info, tol);

[box_pts_1, tmpl_pts_1, tmpl_idxes_1] = abstract_neighbor_spreading_2D(grid_info_1, proxy_info_1);
box_ctr = bin_center(grid_info_1.center_bin, grid_info_1);
disp("test_1: box center: ");
disp(box_ctr);
tmpl_pts_1 = tmpl_pts_1 + box_ctr;  % Shift template points to be centered at the box center

% check the number of interior bins
disp("test_1: grid_info_1.nbin: ");
disp(grid_info_1.nbin);

% Radius in index space
rad = ceil(2 * proxy_info_1.radius / (grid_info_1.nspread * grid_info_1.dx));
disp("test_1: Radius in index space: " + int2str(rad));

% Size check
assert(all(size(box_pts_1) == [2, grid_info_1.nspread^2]), ...
    'box_pts must be [2, nspread^2] for get_grid output');

% Check we are indexing the correct points.
keep_bool = tmpl_idxes_1(1,:) > 0 & tmpl_idxes_1(1,:) <= grid_info_1.ngrid(1) & ...
            tmpl_idxes_1(2,:) > 0 & tmpl_idxes_1(2,:) <= grid_info_1.ngrid(2);

temp_pts = tmpl_pts_1(:, keep_bool);
grid_idxes = (tmpl_idxes_1(1,keep_bool) - 1) * grid_info_1.ngrid(2) + tmpl_idxes_1(2,keep_bool);
grid_pts_idxed = grid_info_1.r(:, grid_idxes);

for i = 1:size(temp_pts, 2)
    dist_x = abs((grid_pts_idxed(1, i) - temp_pts(1, i)));
    dist_x_dx = dist_x / grid_info_1.dx;
    dist_y = abs((grid_pts_idxed(2, i) - temp_pts(2, i)));
    dist_y_dx = dist_y / grid_info_1.dx;
    assert(dist_x < 1e-10, ...
        sprintf('Index %g, computed idx %g, grid point (%g,%g) does not match template point (%g,%g) in x. In units of dx, its %g', ...
                i, grid_idxes(i), grid_pts_idxed(1, i), grid_pts_idxed(2, i), temp_pts(1, i), temp_pts(2, i), dist_x_dx));
    assert(dist_y < 1e-10, ...
        sprintf('Index %g, computed idx %g, grid point (%g,%g) does not match template point (%g,%g) in y. In units of dx, its %g', ...
                i, grid_idxes(i), grid_pts_idxed(1, i), grid_pts_idxed(2, i), temp_pts(1, i), temp_pts(2, i), dist_y_dx));
end

% Figure 1: spreading template in black and box points in red
figure(1);
scatter(tmpl_pts_1(1, :), tmpl_pts_1(2, :), 'kx');
hold on;
scatter(box_pts_1(1, :), box_pts_1(2, :), 'r.');


% Figure 2: source points colored by bin
sort_info = SortInfo(src_info, grid_info_1.dx, grid_info_1.Lbd, grid_info_1.nbin, grid_info_1.nbinpts);
figure(2);
scatter(sort_info.r_srt(1,:), sort_info.r_srt(2,:), 20, sort_info.binid_srt, 'filled');
colormap('parula');
colorbar;
bin0_ctr = bin_center(0, grid_info_1);
hold on;
% Center the spreading template at the center of bin 0 for visualization
scatter(tmpl_pts_1(1, :) + bin0_ctr(1), tmpl_pts_1(2, :) + bin0_ctr(2), 'kx');
% Also add the spreading box
scatter(box_pts_1(1, :) + bin0_ctr(1), box_pts_1(2, :) + bin0_ctr(2), 'ro');
% Also add the proxy ring for bin 0
proxypts = get_ring_points(100, proxy_info_1.radius, bin0_ctr);
plot(proxypts(1,:), proxypts(2,:) , 'k-');
title('Sources colored by bin, with spreading template and box for bin 0');
axis equal;

% Also the neighborhood returned by neighbor_template_2d for bin 0. Plot the 
% nbr_grididxes as text on the figure.
[nbr_binids, nbr_gridpts, nbr_grididxes] = neighbor_template_2d(grid_info_1, proxy_info_1, 0);
nbr_gridpts = nbr_gridpts(:, nbr_grididxes <= grid_info_1.ngrid(1) * grid_info_1.ngrid(2) & nbr_grididxes > 0);
% text(nbr_gridpts(1, :), nbr_gridpts(2, :), arrayfun(@(idx) int2str(idx), nbr_grididxes(nbr_grididxes <= grid_info_1.ngrid(1) * grid_info_1.ngrid(2) & nbr_grididxes > 0), 'UniformOutput', false), 'FontSize', 8);

for i = 1:length(nbr_binids)
    center = bin_center(nbr_binids(i), grid_info_1);
    % Plot the center of the neighboring bin as text on the figure, with the bin id as the text.
    text(center(1), center(2), int2str(nbr_binids(i)), 'FontSize', 8);
end

% scatter(nbr_gridpts(1, :), nbr_gridpts(2, :), 'm.');


% For every bin, count the number of in-bounds points that neighbor_template_2d
% returns. For interior bins this should equal size(tmpl_pts_1, 2).
% n_template = size(tmpl_pts_1, 2);
% ngridpts = grid_info_1.ngrid(1) * grid_info_1.ngrid(2);

% n_bins = grid_info_1.nbin(1) * grid_info_1.nbin(2);
% found_interior = false;
% for bin_idx = 0:(n_bins - 1)
%     [~, ~, nbr_grididxes, ~] = neighbor_template_2d(grid_info_1, proxy_info_1, bin_idx);
%     n_valid = sum(nbr_grididxes <= ngridpts);
%     disp("test_1: For bin_idx " + int2str(bin_idx) + ", nbr_grididxes: " + int2str(size(nbr_grididxes)) + ", n_valid: " + int2str(n_valid) + ", n_template: " + int2str(n_template));
%     if n_valid == n_template
%         found_interior = true;
%         break;
%     end
% end

% assert(found_interior, ...
%     ['No bin found whose valid neighbor grid point count matches ' ...
%      'abstract_neighbor_spreading_2D. Template may be over- or under-sized.']);

%% test_2: 
% Template covers the correct neighbor offsets.
%
% Each point in spreading_template_pts must lie within the spreading box of
% some bin offset (delta_x, delta_y) satisfying the same circular criterion
% used by intersecting_bins_2d.



% dx2    = grid_info.dx;
% nsp2   = grid_info.nspread;
% nbp2   = grid_info.nbinpts;
% rad2   = ceil(2 * proxy_info.radius / (nsp2 * dx2));
% half   = (nsp2 - 1) / 2 * dx2;   % half-width of a spreading box

% for i = 1:size(spreading_template_pts, 2)
%     pt = spreading_template_pts(:, i);
%     % Find which bin offset this point belongs to
%     delta_x = round(pt(1) / (nbp2 * dx2));
%     delta_y = round(pt(2) / (nbp2 * dx2));
%     % Check that (delta_x, delta_y) satisfies the intersection criterion
%     assert(delta_x^2 + delta_y^2 <= rad2^2 + 1e-10, ...
%         sprintf('Template point (%g,%g) belongs to offset (%d,%d) outside neighborhood', ...
%                 pt(1), pt(2), delta_x, delta_y));
%     % Check that the point lies within that bin's spreading box
%     center = [delta_x; delta_y] * nbp2 * dx2;
%     assert(max(abs(pt - center)) <= half + 1e-10, ...
%         sprintf('Template point (%g,%g) not within its bin offset box', pt(1), pt(2)));
% end
