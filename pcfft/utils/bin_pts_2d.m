function [r_sorted, sorted_bin_ids, id_start, sorted_idxes] = bin_pts_2d(r, dx, ngrid, Lbd,  nbin)
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
    disp("bin_pts_2d: Ngrid: ")
    disp(ngrid)

    % if ngrid(1) / nbin is an integer, then we need one extra bin
    if mod(ngrid(1), nbin) == 0
        N_x_bins = ngrid(1)/nbin ;
    else
        N_x_bins = ceil(ngrid(1)/nbin);
    end
    if mod(ngrid(2), nbin) == 0
        N_y_bins = ngrid(2)/nbin + 1;
    else
        N_y_bins = ceil(ngrid(2)/nbin);
    end
    disp("bin_pts_2d: nbin: " + int2str(nbin))
    disp("bin_pts_2d: N_x_bins: " + int2str(N_x_bins))
    disp("bin_pts_2d: N_y_bins: " + int2str(N_y_bins))
    N_bins = N_x_bins * N_y_bins;


    % Find the ID of the bin in the X dim that each 
    % point occupies.
    % NOTE THAT id_x and id_y are zero-indexed!!
    id_x = floor((r(1,:) - Lbd(1) + dx/2) / (nbin * dx));
    % Same for the Y dim.
    id_y = floor((r(2,:) - Lbd(2) + dx/2) / (nbin * dx));
    disp("bin_pts_2d: id_x: " + int2str(id_x))
    disp("bin_pts_2d: id_y: " + int2str(id_y))

    bin_ids = id_x * N_y_bins + id_y;
    disp("bin_pts_2d: bin_ids: " + int2str(bin_ids))
    [sorted_bin_ids, sorted_idxes] = sort(bin_ids);
    % disp(idx_bin_ids)
    % disp(sorted_bin_ids)

    % Sort the points
    r_sorted = r(:, sorted_idxes);

    % Form an array where id_start(i) gives us the index in 
    % r_sorted for the first point with bin idx i.
    % If the bin is empty, id_start(i) = id_start(i-1)

    % We want the slice (id_start(i+1) : id_start(i+2)-1) to give the
    % indices of points in bin i.
    id_start = ones(1,N_bins+1);

    % Loop through sorted_bin_ids and fill in id_start
    current_bin = 0;
    for i = 1:size(r,2)
        bin_i = sorted_bin_ids(i);
        if bin_i > current_bin
            % Fill in all the bins we skipped
            id_start(current_bin+2:bin_i+2) = i;
            current_bin = bin_i;
        end
    end
    % Fill in the rest of the bins
    id_start(current_bin+2:N_bins+1) = size(r, 2) + 1;

