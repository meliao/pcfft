function [box_pts, spreading_template_pts] = abstract_neighbor_spreading_2D(grid_info, proxy_info)
    % Constructs a fictitious spreading box centered at the origin and a
    % spreading template around it, in abstract/relative coordinates.
    %
    % Unlike neighbor_template_2d, this function does not reference any
    % specific bin index or global grid — the result represents the full
    % interior neighborhood for any non-boundary bin.
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
    % box_pts               : array [2, nspread^2]
    %   Coordinates of the spreading box centered at the origin.
    % spreading_template_pts : array [2, n_template_pts]
    %   Coordinates of all unique spreading box points across every
    %   neighboring bin (including the central bin), in origin-relative
    %   coordinates.

    dx      = grid_info.dx;
    nspread = grid_info.nspread;
    nbinpts = grid_info.nbinpts;

    % Spreading box centered at origin: nspread points per dimension.
    box_1d  = dx * (0:nspread-1) - (nspread-1)/2 * dx;
    [X, Y]  = meshgrid(box_1d, box_1d);
    box_pts = [X(:).'; Y(:).'];   % [2, nspread^2]

    % Neighborhood radius in bin-index units — same formula as
    % intersecting_bins_2d (line 28).
    rad = ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));

    % Collect spreading box points for every neighboring bin offset.
    all_pts = zeros(2, 0);
    for delta_x = -rad : rad
        delta_y_max = floor(rad^2 - delta_x^2);
        for delta_y = -delta_y_max : delta_y_max
            shift = [delta_x; delta_y] * nbinpts * dx;
            all_pts = [all_pts, box_pts + shift];
        end
    end

    % Deduplicate while preserving insertion order.
    disp("abstract_neighbor_spreading_2D: Number of points before deduplication: " + int2str(size(all_pts, 2)));
    [~, uid] = unique(all_pts.', 'rows', 'stable');
    spreading_template_pts = all_pts(:, uid);
    disp("abstract_neighbor_spreading_2D: Number of points after deduplication: " + int2str(size(spreading_template_pts, 2)));
end
