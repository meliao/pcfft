function offsets = neighbor_offsets_2d(rad)
    % neighbor_offsets_2d  Integer bin offsets within a circular disk.
    %
    % Returns all integer (delta_x, delta_y) pairs whose Euclidean distance
    % from the origin is at most rad (using ceil rounding on the boundary).
    %
    % Input
    % -----
    % rad : scalar
    %   Neighborhood radius in bin-index units.
    %
    % Output
    % ------
    % offsets : array [2, n]
    %   Each column is one (delta_x; delta_y) offset.

    offsets = zeros(2, 0);
    for delta_x = -rad : rad
        delta_y_max = ceil(sqrt(rad^2 - delta_x^2));
        for delta_y = -delta_y_max : delta_y_max
            offsets = [offsets, [delta_x; delta_y]];
        end
    end
end
