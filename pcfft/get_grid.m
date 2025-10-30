function [grid_info, proxy_info] = get_grid(kernel, src_info, targ_info, ...
        tol, n_nbr)
    % Computes the equispaced spreading grid that will be used for the 
    % precorrected FFT calculation. Takes in a kernel, source points, and 
    % target points, as well as parameters tol and n_nbr. Adaptively finds a 
    % ??? so that <n_nbr> points is approximately the number of source
    % points which are mapped to a single grid point. Then finds grid spacing
    % and discretization of the proxy surface such that accuracy is achieved 
    % by spreading onto a regular grid.
    %
    % Returns information about the regular grid and the proxy points.
    %
    %
    % Inputs kernel : Function mapping (source_pts, target_pts) to a matrix of
    %               size (n_target, n_src).
    %           src_info : struct giving information about the sources
    %               .r : shape (dim, n_src) specifying the location of the
    %                   sources.
    %           targ_info : struct describing the target points.
    %               .r : shape (dim, n_targ) specifying the taget locations.
    %           tol : float specifying absolute error tolerance. Error is 
    %               evaluated at a surface 1.1 * radius of the proxy surface
    %           n_nbr : (optional) int specifying the average number of
    %               interactions that must be done directly. Defaults to
    %               1000.
    %
    % Returns grid_info: struct describing the spreading info with fields
    %               .ngrid : (dim, 1) array containing the number of regular 
    %                       grid points in each dimension.
    %               .Lbd : (dim, 2) array containing the points 
    %                       (xmin, ymin, zmin) and (xmax, ymax, zmax). Describes
    %                       the bounding box of the union of source and target
    %                       points. 
    %               .nspread : integer. Number of spreading points in each 
    %                       dimension necessasry for the requested accuracy.
    %               .rpad : mutliplicative padding factor which extends the 
    %                       size of the bounding box to get rid of artifacts
    %                       when using FFT.
    %               .r : (dim, npts) array containing the points of the padded
    %                       regular grid. In dimension i, there are 
    %                       (rpad * ngrid(i)) + 1 points in the grid.
    %               .dim : integer specifying the dimension of the problem.
    %         proxy_info: struct the proxy points with fields
    %               .n_points_total : integer, the total number of proxy points
    %               .dim : dimension of the problem
    %               .radius : float specifying the distance of proxy
    %                   surface from the center of the spreading box.
    %               .r : (dim, n_points_total) array containing the proxy points

    dim = size(src_info.r, 1);
    % nsrc = size(src_info.r, 2);
    if nargin < 5
        n_nbr = 1000;
    end

    crad = 2;

    % Get the half_sidelen and center of the points to specify the regular grid
    [Lbd, ~] = bounding_box([src_info.r, targ_info.r]);
    halfside = spread_halfside([src_info.r, targ_info.r], n_nbr, crad);

    % get prototype grid for spreading
    [grid_info, proxy_info] = compute_nspread_nproxy(kernel, dim, tol, halfside);

    % get total grid.
    % ngrid is the number of points in the x, y, z direction
    % for the regular grid.
    ngrid = ceil(diff(Lbd, 1, 2) / grid_info.dx);
    rpad = 2;

    if dim == 2
        % Create a regular grid with spacing dx starting at the xmin, ymin point
        % specified by Lbd. Pad this grid 2x the size necessary. The padding
        % is necessary for using the FFT in a later step.
        xx = Lbd(1, 1) + (0:rpad * ngrid(1)) * grid_info.dx;
        yy = Lbd(2, 1) + (0:rpad * ngrid(2)) * grid_info.dx;
        [X, Y] = meshgrid(xx, yy);
        rgrid = [X(:).'; Y(:).'];
    elseif dim == 3
        xx = Lbd(1, 1) + (0:rpad * ngrid(1)) * grid_info.dx;
        yy = Lbd(2, 1) + (0:rpad * ngrid(2)) * grid_info.dx;
        zz = Lbd(3, 1) + (0:rpad * ngrid(3)) * grid_info.dx;
        [X, Y, Z] = meshgrid(xx, yy, zz);
        rgrid = [X(:).'; Y(:).'; Z(:).'];
    end

    % Update grid_info with some more data that we computed.
    grid_info.ngrid = ngrid;
    grid_info.r = rgrid;
    grid_info.Lbd = Lbd;
    grid_info.rpad = rpad;
    grid_info.dim = dim;

end
