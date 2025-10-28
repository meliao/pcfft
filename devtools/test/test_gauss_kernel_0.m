% Ensure things run without error and return the correct size.
addpath('../../pcfft/utils');
n_src = 7;
src_pts = randn(2,n_src);

n_target = 3;
target_pts = randn(2,n_target);

sigma = 0.1;

k = gauss_kernel(src_pts, target_pts,sigma);
% disp(size(k));

assert(all(size(k) == [n_target, n_src]));
assert(all(~isinf(k(:))));
assert(all(~isnan(k(:))));