function k_evals = gauss_kernel(src_pts, target_pts, sigma)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% sigma is a float
% Computes e^{|| src - target||^2 / (2 sigma^2)}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts(1, :) - target_pts(1, :).';
ry = src_pts(2, :) - target_pts(2, :).';

dist_squared = rx.^2 + ry.^2;
k_evals = exp(-dist_squared / (2 * sigma^2));
end





