

n_pts = 100000;
L = 2.0;
Lbd = [-1 -1 1 1];
% r points live on [-1, 1] x [-1, 1]
rng(0);
r = (rand(2, n_pts) - 0.5) * L;
r(2,:) = r(2,:);

% dx = 0.25, so the grid points are at
dx = 0.25;
ngrid = [9 9];
% When we set nbinpts = 3, we expect
% x bins and y bins [-1, -0.25], [-0.25, 0.5], [0.5, 1.]
nbinpts = 3;
nbin = [3 3];

% spoof the GridInfo object. Need nbin, dx, Lbd, nspread, nbinpts, offset, dx 
grid_info = struct;
grid_info.nbin = nbin;
grid_info.dx = dx;
grid_info.Lbd = Lbd;
grid_info.nspread = 2*nbinpts + 1;
grid_info.nbinpts = nbinpts;
pad = ceil((grid_info.nspread - nbinpts)/2);
grid_info.offset = pad * dx - dx/2;
grid_info.dx = dx;

% [a, b, c] = grid_pts_for_bin_2d(0, grid_info);

% Test bin_center function
center_0 = bin_center(0, grid_info);
disp("test_bin_center_0: center of bin 0: ");
disp(center_0);
expected_center_0 = [-0.625; -0.625];
assert(all(abs(center_0 - expected_center_0) < 1e-12));