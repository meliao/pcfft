addpath('../../pcfft/utils');

n = 100;

out = get_sphere_points(n);
assert(all(size(out) == [3 n]));

% Check to ensure that all points are distinct from each other.

diffs_x = out(1, :) - out(1, :).';
diffs_y = out(2, :) - out(2, :).';
diffs_z = out(3, :) - out(3, :).';
diffs = diffs_x.^2 + diffs_y.^2 + diffs_z.^2;
diffs = diffs + eye(size(diffs, 1));
assert(all(diffs(:) > 1e-10));