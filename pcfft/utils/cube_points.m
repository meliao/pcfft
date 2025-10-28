function pts = cube_points(n, half_sidelen, center)
% Returns a discretization of the cube with n points per
% dimension on each face.
if nargin < 3
    center = [0.0 0.0 0.0];
end

xmin = center(1) - half_sidelen;
xmax = center(1) + half_sidelen;
ymin = center(2) - half_sidelen;
ymax = center(2) + half_sidelen;
zmin = center(3) - half_sidelen;
zmax = center(3) + half_sidelen;

% Avoiding the endpoints.
xvals = linspace(xmin, xmax, n);
yvals = linspace(ymin, ymax, n);
zvals = linspace(zmin, zmax, n);

[X, Y] = meshgrid(xvals, yvals);
[~, Z] = meshgrid(yvals, zvals);

X = X(:).';
Y = Y(:).';
Z = Z(:).';

one_vec = ones(1, n^2);

% Face at x = xmin
face_xmin = [xmin * one_vec; Y; Z];
face_xmax = [xmax * one_vec; Y; Z];
face_ymin = [X; ymin * one_vec; Z];
face_ymax = [X; ymax * one_vec; Z];
face_zmin = [X; Y; zmin * one_vec];
face_zmax = [X; Y; zmax * one_vec];

pts = [face_xmin face_xmax face_ymin face_ymax face_zmin face_zmax];
end