function [id_xs, id_ys, binids] = intersecting_bins_2d(bin_idx, grid_info)
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
    % rad = interaction_radius(proxy_info, grid_info);

    % offsets = neighbor_offsets_2d(rad);
    offsets = grid_info.nbr_offsets;
    id_xs = id_x + offsets(1, :);
    id_ys = id_y + offsets(2, :);

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