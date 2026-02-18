function [pts, center, row_idxes] = grid_pts_for_box_3d(bin_idx, grid_info)
    % Imagine a regular grid with bounds [xmin xmax 
    %                                     ymin ymax
    %                                     zmin zmax] = Lbd
    % and grid spacing <dx>. There are [nx ny nz] = ngrid points 
    % in each dimension. There are nspread regular discretization points
    % across the spreading box, and nbinpts regular discretization points for
    % the spreading bin.
    % 
    % NOTE that grid_pts_for_bin_3d is also implemented, make sure you are 
    % using the correct one!
    %
    % This function returns <pts>, the spreading points in the spreading 
    % box. It also returns <center>, the center of the spreading box.
    % It also returns the row indices where <pts> should appear in A_spread.
    %
    % bin_idx = id_x * N_y_bins * N_z_bins + id_y * N_z_bins + id_z;


    dx = grid_info.dx;
    Lbd = grid_info.Lbd;
    ngrid = grid_info.ngrid;

    N_y_bins = grid_info.nbin(2);
    N_z_bins = grid_info.nbin(3);
    nspread = grid_info.nspread;
    nbinpts = grid_info.nbinpts;
    offset = grid_info.offset;

    id_z = mod(bin_idx, N_z_bins);
    id_y = mod(floor(bin_idx / N_z_bins), N_y_bins);
    id_x = floor(bin_idx / (N_y_bins * N_z_bins));

    % Find the grid points corresponding to id_x, id_y, and id_z
    xpts = Lbd(1) - offset + id_x * dx * nbinpts + dx * (0:nspread-1);
    ypts = Lbd(2) - offset + id_y * dx * nbinpts + dx * (0:nspread-1);
    zpts = Lbd(3) - offset + id_z * dx * nbinpts + dx * (0:nspread-1);
    [X,Y,Z] = meshgrid(xpts, ypts, zpts);
    X = permute(X,[3,1,2]);
    Y = permute(Y,[3,1,2]);
    Z = permute(Z,[3,1,2]);
    pts = [X(:).'; Y(:).'; Z(:).'];

    center = [(xpts(1) + xpts(end)) / 2; (ypts(1) + ypts(end)) / 2; (zpts(1) + zpts(end)) / 2];

    % Compute the row indices. 
    % First, we know that the xpts are in position 
    % (id_x - 1) * nbinpts + 1: id_x * nbinpts
    % and the ypts are in position 
    % (id_y - 1) * nbinpts + 1: id_y * nbinpts
    % and the zpts are in position
    % (id_z - 1) * nbinpts + 1: id_z * nbinpts
    % in the padded grid.
    x_positions = id_x * nbinpts + 1: id_x * nbinpts + nspread ;
    y_positions = id_y * nbinpts + 1: id_y * nbinpts + nspread ;
    z_positions = id_z * nbinpts + 1: id_z * nbinpts + nspread ;

    % The padded grid has size ngrid(i)
    % in each dimension i. It loops through the z axis first, then y, then x.
    row_idxes = ones(1, nspread^3);
    for i = 1:nspread
        for j = 1:nspread
            start_row = (i-1) * nspread^2 + (j-1) * nspread + 1;
            end_row = start_row + nspread - 1;
            row_idxes(start_row:end_row) = (x_positions(i)-1) * ngrid(2) * ngrid(3) + (y_positions(j)-1) * ngrid(3) + z_positions ;
        end
    end
end