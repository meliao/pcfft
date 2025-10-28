% Tests return shapes are correct in 2D
addpath("../utils");

n = 17;
rad = 2.0;
pts = get_ring_points(n, rad);

assert(all(size(pts) == [2 n]));