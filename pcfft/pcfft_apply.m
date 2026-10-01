function u = pcfft_apply(mu, A_spread_s, A_spread_t, A_addsub, kern_0hat, sort_info_src, sort_info_targ)
% Compute N-body sum using a precorrected FFT.
%
% Parameters
% ----------
% mu : matrix [opdim(2)*nsrc, 1]
%   Source strengths, in the original source point order.
% A_spread_s : sparse matrix
%   source spreading matrix (see get_spread).
% A_spread_t : sparse matrix
%   target spreading matrix.
% A_addsub : sparse matrix
%   near field corrections matrix (see get_addsub).
% kern_0hat : matrix
%   Fourier transform of background kernel (see get_kernhat).
% sort_info_src : SortInfo
%   sorting of the source points, the second output of get_spread for the
%   sources.
% sort_info_targ : SortInfo
%   sorting of the target points, the second output of get_spread for the
%   targets.
%
% Returns
% -------
% u : matrix
%   potential, in the original target point order.

% The matrices are stored in sorted point order (see get_spread), so move mu
% into that order. Each point has opdim(2) consecutive degrees of freedom.
ps = sort_info_src.ptid_srt(:).';
opdim_s = numel(mu) / numel(ps);
ps = opdim_s*(ps - 1) + (1:opdim_s).';
mu_sorted = mu(ps(:));

sigma_grid = A_spread_s*mu_sorted;
sigma_hat = fftn(reshape(sigma_grid,size(kern_0hat)/2),size(kern_0hat));
u_hat = kern_0hat .* sigma_hat;
ugrid = ifftn(u_hat);
if length(size(kern_0hat)) == 2
    ugrid = ugrid(1:size(kern_0hat,1)/2,1:size(kern_0hat,2)/2);
else
    ugrid = ugrid(1:size(kern_0hat,1)/2,1:size(kern_0hat,2)/2,1:size(kern_0hat,3)/2);
end
u = A_spread_t.'*ugrid(:);
% This version of u is still in the sorted order
u = u + A_addsub*mu_sorted;
% Now put it back into original order
pt = sort_info_targ.ptid_srt(:).';
opdim_t = numel(u) / numel(pt);
pt = opdim_t*(pt - 1) + (1:opdim_t).';
u(pt(:)) = u;
end
