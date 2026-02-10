function [k_evals,grad] = log_kernel(src_pts, target_pts)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% Computes log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts.r(1, :) - target_pts.r(1, :).';
ry = src_pts.r(2, :) - target_pts.r(2, :).';

dist = sqrt(rx.^2 + ry.^2);

k_evals = log(dist);

k_evals(dist < 1e-14) = 0;

if nargout > 1
    grad = zeros(size(rx,1),size(rx,2),2);
    grad(:,:,1) = -rx ./ dist.^2;
    grad(:,:,2) = -ry ./ dist.^2;
    grad(dist < 1e-14,:) = 0;
end
end