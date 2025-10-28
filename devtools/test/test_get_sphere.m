addpath('../../pcfft/utils');

n = 100;

out = get_sphere_points(n);
assert(all(size(out) == [3 n]));