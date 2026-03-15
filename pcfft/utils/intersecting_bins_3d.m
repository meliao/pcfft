function [id_xs, id_ys, id_zs, binids] = intersecting_bins_3d(bin_idx, grid_info, ...
    proxy_info)
    % Given a set of bins which are described by grid_info, and a set of proxy
    % surfaces which are described by proxy_info, return id_xs and id_ys. The 
    % product of these two sets of bins is the set of intersecting bin idxes.
    % 
    % This function may return invalid bin idxes in the first two return values, 
    % i.e. < 0 or >= grid_info.nbin(d). In the third return value, these invalid
    % bin idxes are set to -1.
    % We say that two bins intersect if their proxy surfaces intersect at all.

    N_y_bins = grid_info.nbin(2);
    N_z_bins = grid_info.nbin(3);
    id_z = mod(bin_idx, N_z_bins);
    id_y = mod(floor(bin_idx / N_z_bins), N_y_bins);
    id_x = floor(bin_idx / (N_y_bins * N_z_bins));
    % disp("intersecting_bins_3d: For bin_idx " + int2str(bin_idx) + ...
    %     ", id_x: " + int2str(id_x) + ", id_y: " + int2str(id_y));

    % Which bins are within 2 * radius / (nspread * dx) ?
    id_x_min = id_x - ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
    id_x_max = id_x + ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));

    % id_y_min =  id_y - ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
    % id_y_max =  id_y + ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
    % id_z_min =  id_z - ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
    % id_z_max =  id_z + ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));

    % id_xs = id_x_min:id_x_max;
    % id_ys = id_y_min:id_y_max;
    % id_zs = id_z_min:id_z_max;

    % Radius in index space is 2 * radius / (nspread * dx).
    rad = ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));

    id_xs = [];
    id_ys = [];
    id_zs = [];
    for idx_x = id_x_min:id_x_max
        % Compute max idx_y for this idx_x
        dx_offset = idx_x - id_x;
        idx_y_max = id_y + floor(rad^2 - dx_offset^2 );
        idx_y_min = id_y - floor(rad^2 - dx_offset^2 );

        for idx_y = idx_y_min:idx_y_max
            % Compute max idx_z for this idx_x, idx_y
            dy_offset = idx_y - id_y;
            idx_z_max = id_z + floor(rad^2 - dx_offset^2 - dy_offset^2 );
            idx_z_min = id_z - floor(rad^2 - dx_offset^2 - dy_offset^2 );

            nz = idx_z_max - idx_z_min + 1;
            id_xs = [id_xs, repmat(idx_x, 1, nz)];
            id_ys = [id_ys, repmat(idx_y, 1, nz)];
            id_zs = [id_zs, idx_z_min:idx_z_max];
        end
    end



    % Compute the binids
    binids = id_zs + id_ys * N_z_bins + id_xs * N_y_bins * N_z_bins;

    binids(id_xs<0 | id_xs >= grid_info.nbin(1)) = -1;
    binids(id_ys<0 | id_ys >= grid_info.nbin(2)) = -1;
    binids(id_zs<0 | id_zs >= grid_info.nbin(3)) = -1;
    binids = binids(:);


    binids = binids.';
end