function [box_pts, spreading_template_pts, spreading_template_idxes] = abstract_neighbor_spreading_2D(grid_info, proxy_info)
    % Constructs a spreading box for grid_info.center_bin and a spreading
    % template covering all neighboring bins, in center_bin-relative coordinates.
    %
    % Inputs
    % ------
    % grid_info  : GridInfo
    %   Regular grid and bin structure.
    % proxy_info : ProxyInfo (or equivalent struct with a .radius field)
    %   Proxy geometry used to determine which bins intersect.
    %
    % Outputs
    % -------
    % box_pts                  : array [2, nspread^2]
    %   Coordinates of the spreading box for grid_info.center_bin.
    % spreading_template_pts   : array [2, n_template_pts]
    %   Coordinates of all unique spreading box points across every
    %   neighboring bin (including center_bin), centered at the center of
    %   box_pts.
    % spreading_template_idxes : array [2, n_template_pts]
    %   Grid indices of the spreading template points. Row 1 contains
    %   0-indexed x-grid positions; row 2 contains 0-indexed y-grid positions.
    %   For a target bin, shift these indices by the bin offset to obtain the
    %   corresponding global grid positions.

    % Step 2a: spreading box for center_bin
    [box_pts, box_center] = grid_pts_for_box_2d(grid_info.center_bin, grid_info);

    % Neighborhood radius in bin-index units
    rad = interaction_radius(proxy_info, grid_info);

    % Collect spreading box points for every neighboring bin offset.
    all_pts = zeros(2, 0);
    for delta_x = -rad : rad
        delta_y_max = ceil(sqrt(rad^2 - delta_x^2));
        for delta_y = -delta_y_max : delta_y_max
            shift = [delta_x; delta_y] * grid_info.nbinpts * grid_info.dx;
            all_pts = [all_pts, box_pts + shift];
        end
    end

    % Deduplicate
    % [~, uid] = unique(all_pts.', 'rows');
    % all_unique_pts = all_pts(:, uid);
    all_unique_pts = all_pts;

    % Step 2b: center relative to box_center
    spreading_template_pts = all_unique_pts - box_center;

    % Step 2c: 0-indexed grid positions for each template point
    x_idxes = round((all_unique_pts(1,:) - (grid_info.rmin(1))) / grid_info.dx)+1;
    y_idxes = round((all_unique_pts(2,:) - (grid_info.rmin(2))) / grid_info.dx)+1;
    spreading_template_idxes = [x_idxes; y_idxes];

    % Step 2d: De-duplicate indices
    [~, uid] = unique(spreading_template_idxes.', 'rows');
    spreading_template_pts = spreading_template_pts(:, uid);
    spreading_template_idxes = spreading_template_idxes(:, uid);
end
