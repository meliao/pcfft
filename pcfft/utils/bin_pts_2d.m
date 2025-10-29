function[r_sorted, sorted_bin_ids, id_start] = bin_pts_2d(r, dx, ngrid, Lbd,  nbin)
    % Sort points <r> into a number of bins.
    % Imagine a regular grid with bounds [xmin ymin xmax ymax] = Lbd
    % and grid spacing <dx>. There are [nx ny] = ngrid points 
    % in each dimension.
    % 
    % We want to sort the points into bins which are <nbin> 
    % regular gridpoints across.
    % Create indices for these bins by looping over x first and then y.
    if nargin < 5
        nbin = 1;
    end

    % N_x_bins = ceil(ngrid(1)/nbin);
    N_y_bins = ceil(ngrid(2)/nbin);

    % Find the zero-indexed ID of the bin in the X dim that each 
    % point occupies.
    id_x = floor((r(1,:) - Lbd(1) + dx/2) / (nbin * dx));
    % Same for the Y dim.
    id_y = floor((r(2,:) - Lbd(2) + dx/2) / (nbin * dx)) + 1;

    bin_ids = id_x * N_y_bins + id_y;
    % disp(bin_ids)
    [sorted_bin_ids, sorted_idxes] = sort(bin_ids);
    % disp(idx_bin_ids)
    % disp(sorted_bin_ids)

    % Sort the points
    r_sorted = r(:, sorted_idxes);

    % Form an array where id_start(i) gives us the index in 
    % r_sorted for the first point with bin idx i.
    id_start = zeros(1,sorted_bin_ids(end));
    % Don't forget to fill the first index with 1
    id_start(1) = 1;
    current_idx = 1;
    for i = 1:size(r,2)
        if sorted_bin_ids(i) > current_idx
            id_start(current_idx + 1) = i;
            current_idx = current_idx + 1;
        end
    end
end

