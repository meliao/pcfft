function [u, grad] = planewave(kvec, r)
%PLANEWAVE Evaluate a 2D plane wave and its gradient.
%   [u, grad] = planewave(kvec, r) returns u = exp(i*kvec·r) and
%   grad = grad(u)

rflat = r(:,:);
phase = kvec.' * rflat;
u = exp(1i * phase).';

if nargout > 1
    grad = (1i * u) .* (kvec.');
end
end
