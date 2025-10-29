function [Lbd, center] = bounding_box(pts)
    % Finds a bounding box around pts, and specifies this box via the half_sidelen and 
    % the center of the box.
    % Expects pts to have shape (dim, n_pts)
    % Return value: Lbd is an array of 
    % [xmin xmax
    %  ymin ymax] with shape (2,2)
    %  or [xmin xmax
    %       ymin ymax
    %       zmin zmax] with shape (3,2)
    % Center is the center of the bounding box.

    dim = size(pts, 1);
    
    mins = min(pts,[],2);
    maxs = max(pts,[],2);

    center = (mins + maxs) / 2;
    center = center;

    Lbd = [mins, maxs];


end