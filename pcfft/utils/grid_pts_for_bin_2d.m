function pts = grid_pts_for_bin_2d(bin_idx, dx, ngrid, Lbd, nbin)
    % Imagine a regular grid with bounds [xmin ymin xmax ymax] = Lbd
    % and grid spacing <dx>. There are [nx ny] = ngrid points 
    % in each dimension.
    disp("bin_idx");
    disp(bin_idx);
    % bin_idx = id_x * ngrid(2) + id_y

    N_y_bins = ceil(ngrid(2)/nbin);
    id_y = mod(bin_idx, N_y_bins);
    disp("id_y");
    disp(id_y);
    disp("bin_idx - id_y");
    disp(bin_idx - id_y);
    id_x = floor((bin_idx - id_y)/N_y_bins);

    disp("id_x");
    disp(id_x);

    % Find the grid points corresponding to id_x and id_y
    xpts = Lbd(1) + id_x * dx * nbin + dx * (0:nbin-1);
    disp("xpts")
    disp(xpts);
    ypts = Lbd(2) + id_y * dx * nbin + dx * (0:nbin-1);
    disp(ypts);
    [X,Y] = meshgrid(xpts, ypts);

    pts = [X(:) Y(:)].';

end