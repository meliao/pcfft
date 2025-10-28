function n_per_dim = get_reg_pts_satisfying_fill(src_pts, half_sidelen, fill)
    %

    
    dim = size(src_pts, 1);
    nsrc = size(src_pts, 2);
    % Compute a matrix where m_{ij} = || x_i - x_j||
    rx = src_pts(1, :) - src_pts(1, :).';
    ry = src_pts(2, :) - src_pts(2, :).';
    if dim == 2
        dists = sqrt(rx.^2 + ry.^2);
    else
        rz = src_pts(3, :) - src_pts(3, :).';
        dists = sqrt(rx.^2 + ry.^2 + rz.^2);
    end

    % Now we sort the rows of the dists matrix.
    dists = sort(dists, 1);
    % Drop the first col because that will be ||x_i - x_i|| = 0
    dists = dists(:, 2:nsrc);


    % Find the minimum of column <fill>
    min_dist = min(dists(:, fill));

    % We need 3 dx <= min_dist
    dx = min_dist / 3;

    % Compute number of grid points that will give you the 
    % correct dx
    n_per_dim = (2 * half_sidelen) / dx;

end