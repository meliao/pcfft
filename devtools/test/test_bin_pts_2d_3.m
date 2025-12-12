% Evaluate bin_pts_2d on manually-generated points so we can check the 
% outputs.
addpath(genpath("../../pcfft"));

% Here's the example from the bin_pts_2d documentation.:
% Suppose the points in r live on [-1, 1] x [-0.5, 0.5]
% dx = 0.25, so the grid points are at
% x grid = [-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
% y grid = [-0.5, -0.25, 0.0, 0.25 0.5]
% Then ngrid = [9 5]
% and if we set nbinpts = 3, we expect 
% x bins [-1.125, -0.375], [-0.375, 0.375], [0.375, 1.125]
% y bins [-0.625, 0.125], [0.125, 0.875]

dx = 0.25;
ngrid = [9 5];
Lbd = [-1 1
       -0.5 0.5];
nbinpts = 3;
nbin = [3 2];

% Put points at (-0.9, 0.4) idx_x 0, idx_y 1, bin 1
%               (0.0, 0.0)  idx_x 1, idx_y 0, bin 2
%               (0.8, 0.4)  idx_x 2, idx_y 1, bin 5
r = [-0.9, 0.0, 0.8;
      0.4, 0.0, 0.4];
assert(size(r,1) == 2);
assert(size(r,2) == 3);

[r_sorted, bin_idxes, ~, id_start] = bin_pts_2d(r, dx, Lbd, nbin, nbinpts);

disp("test: r_sorted:")
disp(r_sorted)
disp("test: bin_idxes:")
disp(bin_idxes)
disp("test: id_start:")
disp(id_start)

% We expect bin_idxes = [0, 2, 5]
expected_bin_idxes = [1, 2, 5];
assert(all(bin_idxes == expected_bin_idxes));

% for bin i, we want the slice (id_start(i+1) : id_start(i+2)-1) to give the 
% indices of points in bin i.
% So we expect:
% bin 0: no points -> id_start(1) = 1, id_start(2) = 1
% bin 1: point 1   -> id_start(2) = 1, id_start(3) = 2
% bin 2: point 2   -> id_start(3) = 2, id_start(4) = 3
% bin 3: no points -> id_start(4) = 3, id_start(5) = 3
% bin 4: no points -> id_start(5) = 3, id_start(6) = 3
% bin 5: point 3   -> id_start(6) = 3, id_start(7) = 4
% bin 6: no points -> id_start(7) = 4
expected_id_start = [1, 1, 2, 3, 3, 3, 4];
disp("test: expected_id_start:")
disp(expected_id_start)
assert(all(id_start == expected_id_start));

%% Test case two: new points

% Put points at (-0.9, -0.4) idx_x 0, idx_y 0, bin 0
%               (0.8, 0.4)   idx_x 2, idx_y 1, bin 5
r2 = [-0.9, 0.8;
       -0.4, 0.4];
assert(size(r2,1) == 2);
assert(size(r2,2) == 2);

[r_sorted2, bin_idxes2, ~, id_start2] = bin_pts_2d(r2, dx, Lbd, nbin, nbinpts);

disp("test case 2: r_sorted2:")
disp(r_sorted2)
disp("test case 2: bin_idxes2:")
disp(bin_idxes2)
disp("test case 2: id_start2:")
disp(id_start2)

expected_bin_idxes2 = [0, 5];
assert(all(bin_idxes2 == expected_bin_idxes2));

% For id_start2, we expect:
% bin 0: point 1   -> id_start2(1) = 1, id_start2(2) = 2
% bin 1: no points -> id_start2(2) = 2, id_start2(3) = 2
% bin 2: no points -> id_start2(3) = 2, id_start2(4) = 2
% bin 3: no points -> id_start2(4) = 2, id_start2(5) = 2
% bin 4: no points -> id_start2(5) = 2, id_start2(6) = 2
% bin 5: point 2   -> id_start2(6) = 2, id_start2(7) = 3
% bin 6: no points -> id_start2(7) = 3, id_start2(8) = 3
expected_id_start2 = [1, 2, 2, 2, 2, 2, 3];
disp("test case 2: expected_id_start2:")
disp(expected_id_start2)
assert(all(id_start2 == expected_id_start2));

%% Test case three: new points; now there are 2 points in one bin

% Put points at (-0.9, -0.4) idx_x 0, idx_y 0, bin 0
%               (-0.9, -0.3) idx_x 0, idx_y 0, bin 0
%               (0.0, 0.4)   idx_x 1, idx_y 1, bin 3
r3 = [-0.9, -0.9, 0.0;
       -0.4, -0.3, 0.4];
assert(size(r3,1) == 2);
assert(size(r3,2) == 3);

[r_sorted3, bin_idxes3, ~, id_start3] = bin_pts_2d(r3, dx, Lbd, nbin, nbinpts);
disp("test case 3: r_sorted3:")
disp(r_sorted3)
disp("test case 3: bin_idxes3:")
disp(bin_idxes3)
disp("test case 3: id_start3:")
disp(id_start3)


expected_bin_idxes3 = [0, 0, 3];
assert(all(bin_idxes3 == expected_bin_idxes3));

% For id_start3, we expect:
% bin 0: points 1, 2 -> id_start3(1) = 1, id_start3(2) = 3
% bin 1: no points  -> id_start3(2) = 3, id_start3(3) = 3
% bin 2: no points  -> id_start3(3) = 3, id_start3(4) = 3
% bin 3: no points  -> id_start3(4) = 3, id_start3(5) = 4
% bin 4: no points  -> id_start3(5) = 4, id_start3(6) = 4
% bin 5: point 3   -> id_start3(6) = 4, id_start3(7) = 4
% bin 6: no points  -> id_start3(7) = 4, id_start3(8) = 4
expected_id_start3 = [1, 3, 3, 3, 4, 4, 4];
disp("test case 3: expected_id_start3:")
disp(expected_id_start3)
assert(all(id_start3 == expected_id_start3));

