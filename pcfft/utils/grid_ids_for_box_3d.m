function [row_idxes] = grid_ids_for_box_3d(bin_idx, grid_info)
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
    % This function returns the row indices where <pts> should appear in A_spread.
    %
    % bin_idx = id_x * N_y_bins * N_z_bins + id_y * N_z_bins + id_z;


    ngrid = grid_info.ngrid;

    N_y_bins = grid_info.nbin(2);
    N_z_bins = grid_info.nbin(3);
    nspread = grid_info.nspread;
    nbinpts = grid_info.nbinpts;

    id_z = mod(bin_idx, N_z_bins);
    id_y = mod(floor(bin_idx / N_z_bins), N_y_bins);
    id_x = floor(bin_idx / (N_y_bins * N_z_bins));


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
    row_idxes = (z_positions(:)) + (y_positions(:)-1).'*ngrid(3) + reshape(x_positions-1,1,1,[]) * ngrid(2) * ngrid(3);
    row_idxes = row_idxes(:).';
end