function u = pcfft_apply(sigma, A_spread_s, A_spread_t, A_addsub, kern_0hat)
% Compute N-body sum using a precorrected FFT.
%
% Parameters
% ----------
% sigma : matrix
%   TODO
% A_spread_s : sparse matrix
%   source spreading matrix (see get_spread).
% A_spread_t : sparse matrix
%   target spreading matrix.
% A_addsub : sparse matrix
%   near field corrections matrix (see get_addsub).
% kern_0hat : matrix
%   Fourier transform of background kernel (see get_kernhat).
%
% Returns
% -------
% u : matrix
%   potential.


sigma_grid = A_spread_s*sigma;
sigma_hat = fftn(reshape(sigma_grid,size(kern_0hat)/2),size(kern_0hat));
u_hat = kern_0hat .* sigma_hat;
ugrid = ifftn(u_hat);
if length(size(kern_0hat)) == 2
    ugrid = ugrid(1:size(kern_0hat,2)/2,1:size(kern_0hat,1)/2);
else
    ugrid = ugrid(1:size(kern_0hat,1)/2,1:size(kern_0hat,2)/2,1:size(kern_0hat,3)/2);
end
u = A_spread_t.'*ugrid(:) + A_addsub*sigma;

end
