function k_evals = log_kernel(src_pts, target_pts)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% Computes log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts(1, :) - target_pts(1, :).';
ry = src_pts(2, :) - target_pts(2, :).';

dist = sqrt(rx.^2 + ry.^2);

k_evals = log(dist);

k_evals(dist < 1e-14) = 0;
end