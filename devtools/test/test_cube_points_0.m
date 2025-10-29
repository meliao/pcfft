addpath('../../pcfft/utils');

n = 3;
half_sidelen = 0.5;

pts = cube_points(n, half_sidelen);

assert(all(size(pts) == [3 6*n^2]))