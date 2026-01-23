function k_evals = one_over_r_kernel2D(src_pts, target_pts)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% sigma is a float
% Computes 1  / || src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts.r(1, :) - target_pts.r(1, :).';
ry = src_pts.r(2, :) - target_pts.r(2, :).';

dist = sqrt(rx.^2 + ry.^2);
k_evals = 1 ./ dist;
k_evals(dist < 1e-14) = 0;
end