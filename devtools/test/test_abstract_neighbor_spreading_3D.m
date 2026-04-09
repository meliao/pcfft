% Tests for abstract_neighbor_spreading_2D.
addpath(genpath('../../pcfft'));

%% test_0: 
% Plot a basic spreading template to visually check it.

dim = 3;
Lbd = [-1 1; -1 1; -1 1];

nsrc = 1000;
ntarg = 1000;
rng(0);
src_info = struct;
src_info.r = (rand(dim, nsrc) - 0.5) * 2;
targ_info = struct;
targ_info.r = (rand(dim, ntarg) - 0.5) * 2;

[grid_info, proxy_info] = get_grid(@one_over_r_kernel, src_info, targ_info, 1e-6);

[box_pts, spreading_template_pts, spreading_template_idxes] = abstract_neighbor_spreading_3D(grid_info, proxy_info);

% Scatter plot spreading_template_pts
scatter3(spreading_template_pts(1, :), spreading_template_pts(2, :), spreading_template_pts(3, :), 'filled');

% Also scatter the gridpoints in a different color
hold on;
scatter3(grid_info.r(1, :), grid_info.r(2, :), grid_info.r(3, :), 'k.');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Spreading Template Points');

%% test_1: 
% Basic size and centering checks on a hand-constructed grid.
%
% dx=0.5, nbinpts=2, nspread=4, so each spreading box is 4x4=16 pts.

dim = 3;
Lbd = [-1 1; -1 1; -1 1];
dx = 0.5;
nbinpts = 2;
nspread = 4;
grid_info = GridInfo(Lbd, dx, nspread, nbinpts, dim, -1);

proxy_info = struct;
proxy_info.radius = 0.5;

[box_pts, spreading_template_pts, spreading_template_idxes] = abstract_neighbor_spreading_3D(grid_info, proxy_info);




% box_pts must be [3, nspread^3]
assert(all(size(box_pts) == [3, nspread^3]), ...
    'box_pts must be [3, nspread^3]');

% box_pts must be centered at center_bin's spreading box center
[~, expected_box_center] = grid_pts_for_box_3d(grid_info.center_bin, grid_info);
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

% spreading_template_pts must have at least nspread^3 points
assert(size(spreading_template_pts, 2) >= nspread^3, ...
    'spreading_template_pts must have at least nspread^3 points');

% spreading_template_pts must be separated by about dx
dists = sqrt(sum(diff(spreading_template_pts, 1, 2).^2, 1));
min_dist = min(dists);
assert(min_dist + 1e-12 >= dx, ...
    sprintf('Minimum distance between spreading_template_pts is %g, less than dx=%g', min_dist, dx));



% All indices must be unique
[~, uid] = unique(spreading_template_idxes.', 'rows');
assert(length(uid) == size(spreading_template_idxes, 2), ...
    'spreading_template_idxes must have no duplicate indices');