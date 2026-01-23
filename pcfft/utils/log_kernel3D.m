function k_evals = log_kernel3D(src_pts, target_pts)
% src_pts has shape (3, M)
% target_pts has shape (3, N)
% Computes log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts.r(1, :) - target_pts.r(1, :).';
ry = src_pts.r(2, :) - target_pts.r(2, :).';
rz = src_pts.r(3, :) - target_pts.r(3, :).';

dist = sqrt(rx.^2 + ry.^2 + rz.^2);

k_evals = log(dist);

k_evals(dist < 1e-14) = 0;
end