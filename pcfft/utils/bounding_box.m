function [Lbd, center] = bounding_box(pts)
    % Finds a bounding box around pts, and specifies this box via the half_sidelen and 
    % the center of the box.
    % Expects pts to have shape (dim, n_pts)

    dim = size(pts, 1);
    
    mins = min(pts.');
    maxs = max(pts.');

    center = (mins + maxs) / 2;
    center = center.';

    Lbd = [mins, maxs];


end