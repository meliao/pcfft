addpath('../../pcfft/utils');

% Bounding box is
% (0.0, 0.1), (0.1, 0.1) (0.1, 0.2) (0.0, 0.2)
% so half_sidelen = 0.05
% and center = (0.05, 0.15)
pts = [0.0 0.1 0.05
        0.1 0.2 0.1];
% disp(size(pts));

[half_sidelen, center] = bounding_box(pts);

expected_half_sidelen = 0.05;
% disp(half_sidelen);
assert(half_sidelen == expected_half_sidelen);
expected_center = [0.05 0.15].';
% disp(center);
% disp(expected_center);
assert(all(size(center) == [2 1]));
diffs = expected_center(:) -center(:);
assert(all(diffs < 1e-15));


% 3D example
pts = [0.0 1.0 0.0 0.0
        1.0 0.0 0.0 0.0
        0.5 0.0 0.0 0.1];
[half_sidelen, center] = bounding_box(pts);

expected_half_sidelen = 0.5;
expected_center = [0.5 0.5 0.25].';
% disp(half_sidelen);
assert(half_sidelen == expected_half_sidelen);
assert(all(size(center) == [3 1]));
diffs = expected_center(:) -center(:);
assert(all(diffs < 1e-15));