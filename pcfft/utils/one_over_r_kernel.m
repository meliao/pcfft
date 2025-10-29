function k_evals = one_over_r_kernel(src_pts, target_pts)
% src_pts has shape (3, M)
% target_pts has shape (3, N)
% sigma is a float
% Computes 1  / || src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts(1, :) - target_pts(1, :).';
ry = src_pts(2, :) - target_pts(2, :).';
rz = src_pts(3, :) - target_pts(3, :).';

dist = sqrt(rx.^2 + ry.^2 + rz.^2);
k_evals = 1 ./ dist;
k_evals(dist < 1e-14) = 0;
end