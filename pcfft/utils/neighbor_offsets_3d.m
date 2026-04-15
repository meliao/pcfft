function offsets = neighbor_offsets_3d(rad)
    % neighbor_offsets_3d  Integer bin offsets within a spherical ball.
    %
    % Returns all integer (delta_x, delta_y, delta_z) triples whose Euclidean
    % distance from the origin is at most rad (using ceil rounding on the boundary).
    %
    % Input
    % -----
    % rad : scalar
    %   Neighborhood radius in bin-index units.
    %
    % Output
    % ------
    % offsets : array [3, n]
    %   Each column is one (delta_x; delta_y; delta_z) offset.

    offsets = zeros(3, 0);
    for delta_x = -rad : rad
        delta_y_max = floor(sqrt(rad^2 - delta_x^2));
        for delta_y = -delta_y_max : delta_y_max
            r2 = rad^2 - delta_x^2 - delta_y^2;
            if r2 < 0; continue; end
            delta_z_max = floor(sqrt(r2));
            delta_zs = -delta_z_max : delta_z_max;
            offsets = [offsets, [delta_x * ones(1, numel(delta_zs)); delta_y * ones(1, numel(delta_zs)); delta_zs]];
        end
    end
end
