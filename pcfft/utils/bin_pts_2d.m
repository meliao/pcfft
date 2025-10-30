function [r_sorted, sorted_bin_ids, id_start] = bin_pts_2d(r, dx, ngrid, Lbd,  nbin)
    % Sort points <r> into a number of bins.
    % Imagine a regular grid with bounds [xmin ymin xmax ymax] = Lbd
    % and grid spacing <dx>. There are [nx ny] = ngrid points 
    % in each dimension.
    % 
    % We want to sort the points into bins which are <nbin> 
    % regular gridpoints across.
    % Create indices for these bins by looping over x first and then y.
    %
    %
    % Example:
    % Suppose the points in r live on [-1, 1] x [-0.5, 0.5]
    % dx = 0.25, so the grid points are at
    % x grid = [-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
    % y grid = [-0.5, -0.25, 0.0, 0.25 0.5]
    % Then ngrid = [9 5]
    % and if we set nbin = 3, we expect 
    % x bins [-1.125, -0.375], [-0.375, 0.375], [0.375, 1.125]
    % y bins [-0.625, 0.125], [0.125, 0.875]
    if nargin < 5
        nbin = 1;
    end

    N_x_bins = ceil(ngrid(1)/nbin);
    N_y_bins = ceil(ngrid(2)/nbin);
    disp("N_x_bins: " + int2str(N_x_bins))
    disp("N_y_bins: " + int2str(N_y_bins))
    N_bins = N_x_bins * N_y_bins;


    % Find the ID of the bin in the X dim that each 
    % point occupies.
    % NOTE THAT id_x and id_y are zero-indexed!!
    id_x = floor((r(1,:) - Lbd(1) + dx/2) / (nbin * dx));
    % Same for the Y dim.
    id_y = floor((r(2,:) - Lbd(2) + dx/2) / (nbin * dx));

    bin_ids = id_x * N_y_bins + id_y;
    % disp(bin_ids)
    [sorted_bin_ids, sorted_idxes] = sort(bin_ids);
    % disp(idx_bin_ids)
    % disp(sorted_bin_ids)

    % Sort the points
    r_sorted = r(:, sorted_idxes);

    % Form an array where id_start(i) gives us the index in 
    % r_sorted for the first point with bin idx i.
    % If the bin is empty, id_start(i) = id_start(i-1)
    id_start = zeros(1,N_bins+1);
    ibin = 1;
    id_start(1) = 1;
    id_start(end) = size(r, 2) + 1;

    % disp("sorted_bin_ids");
    % disp(sorted_bin_ids);


    for i = 1:size(r,2)
        if sorted_bin_ids(i) > ibin
            ibinold = ibin;
            ibin = sorted_bin_ids(i);
            id_start(ibinold+1:ibin) = i;
        end
    end
end
