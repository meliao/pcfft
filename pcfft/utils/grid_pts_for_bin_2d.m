function [pts, center, row_idxes] = grid_pts_for_bin_2d(bin_idx, grid_info, nbin)
    % Imagine a regular grid with bounds [xmin ymin xmax ymax] = Lbd
    % and grid spacing <dx>. There are [nx ny] = ngrid points 
    % in each dimension.
    %
    % Also returns the row indices where <pts> should appear in A_spread.
    %
    % bin_idx = id_x * N_y_bins + id_y

    dx = grid_info.dx;
    Lbd = grid_info.Lbd;
    rpad = grid_info.rpad;
    ngrid = grid_info.ngrid;

    if mod(ngrid(2), nbin) == 0
        N_y_bins = ngrid(2)/nbin + 1;
    else
        N_y_bins = ceil(ngrid(2)/nbin);
    end
    id_y = mod(bin_idx, N_y_bins);
    id_x = floor((bin_idx - id_y)/N_y_bins);


    % Find the grid points corresponding to id_x and id_y
    xpts = Lbd(1) + id_x * dx * nbin + dx * (0:nbin-1);
    ypts = Lbd(2) + id_y * dx * nbin + dx * (0:nbin-1);
    [X,Y] = meshgrid(xpts, ypts);

    pts = [X(:) Y(:)].';

    center = [(xpts(1) + xpts(end)) / 2; (ypts(1) + ypts(end)) / 2];


    % Compute the row indices
    % First, we know that the xpts are in position 
    % (id_x - 1) * nbin + 1: id_x * nbin
    % and the ypts are in position 
    % (id_y - 1) * nbin + 1: id_y * nbin
    % in the padded grid.
    x_positions = id_x * nbin + 1: (id_x + 1) * nbin;
    y_positions = id_y * nbin + 1: (id_y + 1) * nbin;


    % The padded grid has size rpad * ngrid(i)
    % in each dimension i. It loops through the y axis first.
    % Remember there are (rpad * ngrid(2) + 1) pts in the y direction
    row_idxes = ones(1, nbin^2);
    for i = 1:nbin
        start_row = (i-1) * nbin + 1;
        end_row = i * nbin;

        row_idxes(start_row:end_row) = (x_positions(i)-1) * (rpad * ngrid(2)+1) + y_positions ;
    end

end