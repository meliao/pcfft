function u = planewave(kvec, pts)
% Evaluate a plane wave exp(i * kvec @ pts).
%
% Parameters
% ----------
% kvec : matrix
%   TODO
% pts : matrix
%   TODO
%
% Returns
% -------
% u : matrix
%   TODO

kvec = kvec(:);
u = exp(1i * (kvec.' * pts));

end
