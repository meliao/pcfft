function bin_idxes = intersecting_bins_2d(bin_idx, grid_info, proxy_info)
    % Given a set of bins which are described by grid_info, and a set of proxy
    % surfaces which are described by proxy_info, return a list of bins which
    % intersect with the bin which has idx = bin_idx.
    %
    % We say that two bins intersect if their proxy surfaces intersect at all.

    center_this_bin = bin_center(bin_idx, grid_info);
    rad = proxy_info.radius;

    N_bins = grid_info.nbin(1) * grid_info.nbin(2);

    % Loop through the other bins, find their centers, and evaluate whether 
    % the centers are within 2 * rad from each other.
    bin_idxes = [];
    for i = 0:N_bins
        center_i = bin_center(i, grid_info);
        dist = norm(center_this_bin - center_i);
        if dist <= 2 * rad
            bin_idxes = [bin_idxes, i];
        end
    end
end