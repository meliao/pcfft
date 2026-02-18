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

    id_y_min =  id_y - ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
    id_y_max =  id_y + ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
    id_z_min =  id_z - ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));
    id_z_max =  id_z + ceil(2 * proxy_info.radius / (grid_info.nspread * grid_info.dx));

    id_xs = id_x_min:id_x_max;
    id_ys = id_y_min:id_y_max;
    id_zs = id_z_min:id_z_max;

    % Compute the binids
    binids = zeros(length(id_xs) * length(id_ys) * length(id_zs), 1 );
    for i = 1:length(id_xs)
        for j = 1:length(id_ys)
            for k = 1:length(id_zs)
                current_idx = (i-1) * length(id_ys) * length(id_zs) + (j-1) * length(id_zs) + k;

                % If it's an invalid binid, set it to -1
                if id_xs(i) < 0 || id_xs(i) >= grid_info.nbin(1) || ...
                id_ys(j) < 0 || id_ys(j) >= grid_info.nbin(2) || ...
                id_zs(k) < 0 || id_zs(k) >= grid_info.nbin(3)
                    binids(current_idx) = -1;
                else
                    binids(current_idx) = ...
                    id_xs(i) * N_y_bins * N_z_bins + id_ys(j) * N_z_bins + id_zs(k);
                end
            end
        end
    end
    binids = binids.';
end