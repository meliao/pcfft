% Tests return shapes are correct in 2D and 3D
addpath('../../pcfft/utils');


dim = 2;
n = 13;
half_sidelen = 0.5;
out = get_regular_grid(n, half_sidelen, dim);
assert(all(size(out) == [dim, n^dim]));
assert(all(abs(out(:)) <= 1.0));