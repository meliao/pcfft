function [pts, center, row_idxes] = grid_pts_for_box_2d(bin_idx, grid_info)
    % Imagine a regular grid with bounds [xmin ymin xmax ymax] = Lbd
    % and grid spacing <dx>. There are [nx ny] = ngrid points 
    % in each dimension. There are nspread regular discretization points
    % across the spreading box, and nbinpts regular discretization points for
    % the spreading bin.
    % 
    %
    % This function returns <pts>, the spreading points in the spreading 
    % box. It also returns <center>, the center of the spreading box.
    % It also returns the row indices where <pts> should appear in A_spread.
    %
    % bin_idx = id_x * N_y_bins + id_y

    dx = grid_info.dx;
    Lbd = grid_info.Lbd;
    ngrid = grid_info.ngrid;

    N_x_bins = grid_info.nbin(1);
    N_y_bins = grid_info.nbin(2);
    nspread = grid_info.nspread;
    nbinpts = grid_info.nbinpts;
    offset = grid_info.offset;


    id_y = mod(bin_idx, N_y_bins);
    id_x = floor((bin_idx - id_y)/N_y_bins);
    % disp("grid_pts_for_box_2d: id_y: ")
    % disp(id_y)
    
    % disp("grid_pts_for_box_2d: id_x: ")
    % disp(id_x)
    

    % Find the grid points corresponding to id_x and id_y
    xpts = Lbd(1) - offset + id_x * dx * nbinpts + dx * (0:nspread-1);
    ypts = Lbd(2) - offset + id_y * dx * nbinpts + dx * (0:nspread-1);
    [X,Y] = meshgrid(xpts, ypts);

    pts = [X(:) Y(:)].';

    center = [(xpts(1) + xpts(end)) / 2; (ypts(1) + ypts(end)) / 2];
    % disp("grid_pts_for_box_2d: center: ");

    % disp(center);

    % Compute the row indices
    % First, we know that the xpts are in position 
    % (id_x - 1) * nbinpts + 1: id_x * nbinpts
    % and the ypts are in position 
    % (id_y - 1) * nbinpts + 1: id_y * nbinpts
    % in the padded grid.
    x_positions = id_x * nbinpts + 1: id_x * nbinpts + nspread ;
    y_positions = id_y * nbinpts + 1: id_y * nbinpts + nspread ;

    % The padded grid has size ngrid(i)
    % in each dimension i. It loops through the y axis first.
    % Remember there are ngrid(2) pts in the y direction
    row_idxes = ones(1, nspread^2);
    for i = 1:nspread
        start_row = (i-1) * nspread + 1;
        end_row = start_row + nspread - 1;

        row_idxes(start_row:end_row) = (x_positions(i)-1) * ngrid(2) + y_positions ;
    end

end