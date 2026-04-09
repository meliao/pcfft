function [box_pts, spreading_template_pts, spreading_template_idxes] = abstract_neighbor_spreading_3D(grid_info, proxy_info)
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
    % box_pts                  : array [3, nspread^2]
    %   Coordinates of the spreading box for grid_info.center_bin.
    % spreading_template_pts   : array [3, n_template_pts]
    %   Coordinates of all unique spreading box points across every
    %   neighboring bin (including center_bin), centered at the center of
    %   box_pts.
    % spreading_template_idxes : array [3, n_template_pts]
    %   Grid indices of the spreading template points. Row 1 contains
    %   0-indexed x-grid positions; row 2 contains 0-indexed y-grid positions.
    %   Row 3 contains 0-indexed z-grid positions. For a target bin, shift
    %   these indices by the bin offset to obtain the corresponding global
    %   grid positions.


    % First, get the spreading box for the center bin.
    [box_pts, box_center] = grid_pts_for_box_3d(grid_info.center_bin, grid_info);

    % Neighbor radius in bin-index units
    rad = interaction_radius(proxy_info, grid_info);

    % Collect spreading box points for every neighboring bin offset.
    all_pts = zeros(3, 0);
    offsets = neighbor_offsets_3d(rad);
    for k = 1 : size(offsets, 2)
        shift = offsets(:, k) * grid_info.nbinpts * grid_info.dx;
        all_pts = [all_pts, box_pts + shift];
    end

    % Center relative to box_center
    spreading_template_pts = all_pts - box_center;

    % Get 0-indexed grid positions for each template points
    x_idxes = round((all_pts(1,:) - (grid_info.rmin(1))) / grid_info.dx)+1;
    y_idxes = round((all_pts(2,:) - (grid_info.rmin(2))) / grid_info.dx)+1;
    z_idxes = round((all_pts(3,:) - (grid_info.rmin(3))) / grid_info.dx)+1;
    spreading_template_idxes = [x_idxes; y_idxes; z_idxes];

    % Deduplicate indices
    [spreading_template_idxes, unique_idx] = unique(spreading_template_idxes.', 'rows');
    spreading_template_pts = spreading_template_pts(:, unique_idx);

end