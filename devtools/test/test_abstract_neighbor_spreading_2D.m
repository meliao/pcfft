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

% box_pts must be centered at the origin
box_center = mean(box_pts, 2);
assert(norm(box_center) < 1e-12, ...
    'box_pts must be centered at the origin');

% No duplicate points in spreading_template_pts
[~, uid] = unique(spreading_template_pts.', 'rows');
assert(length(uid) == size(spreading_template_pts, 2), ...
    'spreading_template_pts must have no duplicate points');

% box_pts must be a subset of spreading_template_pts (central bin is included)
for i = 1:size(box_pts, 2)
    dists = sqrt(sum((spreading_template_pts - box_pts(:, i)).^2, 1));
    assert(min(dists) < 1e-12, ...
        sprintf('box_pts column %d not found in spreading_template_pts', i));
end

% spreading_template_pts must have at least nspread^2 points
assert(size(spreading_template_pts, 2) >= nspread^2, ...
    'spreading_template_pts must have at least nspread^2 points');

% Spacing between adjacent points must be dx everywhere
% (all points lie on the same dx-spaced grid as box_pts)
% x_vals = unique(round(spreading_template_pts(1, :) / dx));
% y_vals = unique(round(spreading_template_pts(2, :) / dx));
% x_spacings = diff(sort(x_vals));
% y_spacings = diff(sort(y_vals));
% disp("x_vals: ");
% disp(x_vals);
% disp("y_vals: ");
% disp(y_vals);
% disp("x_spacings: ");
% disp(x_spacings);
% disp("y_spacings: ");
% disp(y_spacings);
% assert(all(abs(x_spacings - 1) < 1e-10), ...
%     'x-coordinates of spreading_template_pts must be evenly spaced by dx');
% assert(all(abs(y_spacings - 1) < 1e-10), ...
%     'y-coordinates of spreading_template_pts must be evenly spaced by dx');

%% test_1: 
% Template size matches an interior bin from neighbor_template_2d.
%
% For a grid large enough that interior bins exist, the number of valid
% (in-bounds) points returned by neighbor_template_2d for an interior bin
% must equal size(spreading_template_pts, 2).

rng(0);
n_src = 200;
n_targ = 150;
tol = 1e-8;

scale = 4.0;

src_info = struct;
src_info.r = (rand(2, n_src) - 0.5) * scale;
src_info.weights = rand(n_src, 1);

targ_info = struct;
targ_info.r = (rand(2, n_targ) - 0.5) * scale;

[grid_info_1, proxy_info_1] = get_grid(@log_kernel, src_info, targ_info, tol, 50);

[box_pts_1, tmpl_pts_1] = abstract_neighbor_spreading_2D(grid_info_1, proxy_info_1);

% check the number of interior bins
disp("test_1: grid_info_1.nbin: ");
disp(grid_info_1.nbin);

% Radius in index space
rad = ceil(2 * proxy_info_1.radius / (grid_info_1.nspread * grid_info_1.dx));

disp("test_1: Radius in index space: " + int2str(rad));
% Size check
assert(all(size(box_pts_1) == [2, grid_info_1.nspread^2]), ...
    'box_pts must be [2, nspread^2] for get_grid output');

% For every bin, count the number of in-bounds points that neighbor_template_2d
% returns. For interior bins this should equal size(tmpl_pts_1, 2).
n_template = size(tmpl_pts_1, 2);
ngridpts = grid_info_1.ngrid(1) * grid_info_1.ngrid(2);

n_bins = grid_info_1.nbin(1) * grid_info_1.nbin(2);
found_interior = false;
for bin_idx = 0:(n_bins - 1)
    [~, ~, nbr_grididxes, ~] = neighbor_template_2d(grid_info_1, proxy_info_1, bin_idx);
    n_valid = sum(nbr_grididxes <= ngridpts);
    % disp("test_1: For bin_idx " + int2str(bin_idx) + ", n_valid: " + int2str(n_valid) + ", n_template: " + int2str(n_template));
    if n_valid == n_template
        found_interior = true;
        break;
    end
end

assert(found_interior, ...
    ['No bin found whose valid neighbor grid point count matches ' ...
     'abstract_neighbor_spreading_2D. Template may be over- or under-sized.']);

%% test_2: 
% Template covers the correct neighbor offsets.
%
% Each point in spreading_template_pts must lie within the spreading box of
% some bin offset (delta_x, delta_y) satisfying the same circular criterion
% used by intersecting_bins_2d.

dx2    = grid_info.dx;
nsp2   = grid_info.nspread;
nbp2   = grid_info.nbinpts;
rad2   = ceil(2 * proxy_info.radius / (nsp2 * dx2));
half   = (nsp2 - 1) / 2 * dx2;   % half-width of a spreading box

for i = 1:size(spreading_template_pts, 2)
    pt = spreading_template_pts(:, i);
    % Find which bin offset this point belongs to
    delta_x = round(pt(1) / (nbp2 * dx2));
    delta_y = round(pt(2) / (nbp2 * dx2));
    % Check that (delta_x, delta_y) satisfies the intersection criterion
    assert(delta_x^2 + delta_y^2 <= rad2^2 + 1e-10, ...
        sprintf('Template point (%g,%g) belongs to offset (%d,%d) outside neighborhood', ...
                pt(1), pt(2), delta_x, delta_y));
    % Check that the point lies within that bin's spreading box
    center = [delta_x; delta_y] * nbp2 * dx2;
    assert(max(abs(pt - center)) <= half + 1e-10, ...
        sprintf('Template point (%g,%g) not within its bin offset box', pt(1), pt(2)));
end
