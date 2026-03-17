function [id_xs, id_ys, binids] = intersecting_bins_2d(bin_idx, grid_info, ...
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
    id_y = mod(bin_idx, N_y_bins);
    id_x = floor((bin_idx - id_y)/N_y_bins);
    % disp("intersecting_bins_2d: For bin_idx " + int2str(bin_idx) + ...
    %     ", id_x: " + int2str(id_x) + ", id_y: " + int2str(id_y));

    % Radius in index space 
    rad = ceil(2 * proxy_info.radius / (grid_info.nbinpts * grid_info.dx));

    % Which bins are within 2 * radius / (nspread * dx) ?
    id_x_min = id_x - rad;
    id_x_max = id_x + rad;


    id_xs = [];
    id_ys = [];
    % Loop over idx_x
    for idx_x = id_x_min:id_x_max
        
        % Compute max idx_y for this idx_x
        dx_offset = idx_x - id_x;
        idx_y_max = id_y + ceil(sqrt(rad^2 - dx_offset^2));
        idx_y_min = id_y - ceil(sqrt(rad^2 - dx_offset^2));

        ny = idx_y_max - idx_y_min + 1;
        id_xs = [id_xs, repmat(idx_x, 1, ny)];
        id_ys = [id_ys, idx_y_min:idx_y_max];
    end

    % Compute the binids
    binids = id_ys(:) + id_xs(:) * N_y_bins;
    binids(id_xs < 0 | id_xs >= grid_info.nbin(1) | ...
           id_ys < 0 | id_ys >= grid_info.nbin(2)) = -1;
    % for i = 1:length(id_xs)
    %     % If it's an invalid binid, set it to -1
    %     if id_xs(i) < 0 || id_xs(i) >= grid_info.nbin(1) || ...
    %        id_ys(i) < 0 || id_ys(i) >= grid_info.nbin(2)
    %            binids(i) = -1;
    %         else
    %             binids(i) = ...
    %             id_xs(i) * N_y_bins + id_ys(i);
    %         end
    % end
    binids = binids.';
end