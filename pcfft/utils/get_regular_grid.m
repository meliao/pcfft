function pts = get_regular_grid(n, half_sidelen, dim, center)
% n is an integer specifying number of regular points per dimension.
% half_sidelen is a float.
% dim is an integer specifying dimension of the regular grid.
% center is optional to specify center of the cell
% 
% output has shape (dim, n^dim)
if nargin < 4
    center = zeros(dim, 1);
end


xvals = linspace(center(1) - half_sidelen, center(1) + half_sidelen, n);
yvals = linspace(center(2) - half_sidelen, center(2) + half_sidelen, n);
if dim == 2
    [X, Y] = meshgrid(xvals, yvals);
    pts = [X(:).'; Y(:).'];
else
    zvals = linspace(center(3) - half_sidelen, center(3) + half_sidelen, n);
    [X, Y, Z] = meshgrid(xvals, yvals, zvals);
    pts = [X(:).'; Y(:).'; Z(:).'];
end