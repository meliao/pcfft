function rad = interaction_radius(proxy_info, grid_info)
    % Returns the interaction radius in bin-index units, which is used to
    % determine which neighboring bins to spread to.
    %
    % proxy_info : struct with fields
    %   radius : scalar
    %     The radius of the proxy points in physical units.
    %
    % grid_info : struct with fields
    %   nspread : scalar
    %     The number of points in the spreading box along one dimension.
    %   dx : scalar
    %     The grid spacing in physical units.
    
    rad = ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
end