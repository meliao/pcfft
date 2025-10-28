addpath('../../pcfft/utils');

% 2D case
n = 3;
half_sidelen = 0.5;
n_src = 200;
src_pts = (rand(2, n_src)-0.5) * half_sidelen;
fill = 10;

n_per_dim = get_reg_pts_satisfying_fill(src_pts, half_sidelen, fill);


% 3D case
half_sidelen = 0.5;
n_src = 200;
src_pts = (rand(3, n_src)-0.5) * half_sidelen;
fill = 10;

n_per_dim = get_reg_pts_satisfying_fill(src_pts, half_sidelen, fill);
