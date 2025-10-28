function [spread_info, proxy_info, rgrid] = get_grid(kernel, src_info, targ_info, tol, n_nbr)
    % Computes the number of regular gridpoints on a box or cube with size
    % <half_sidelen> needed to approximate an evaluation of the kernel at a
    % point at least <radius> away.
    %
    %
    % Inputs kernel : Function mapping (source_pts, target_pts) to a matrix of
    %               size (n_target, n_src).
    %           src_info : struct giving information about the sources
    %               .r : shape (dim, n_src) specifying the location of the 
    %                   sources.
    %               .weights : shape (n_src, 1) specifying the charges.
    %           targ_info : struct describing the target points.
    %               .radius : minimum distance from the center of the bounding box of 
    %                   source points to the target points.
    %           tol : float specifying absolute error tolerance. Evaluated at a surface
    %               which is 1.1 * radius of the proxy surface
    %           n_nbr : (optional) int specifying the average number of
    %               interactions that must be done directly. Defaults to
    %               1000
    %
    % Returns grid_info: struct with fields
    %               .n_per_dim : integer describing the number of points per dimension
    %               .dim : integer specifying the dimension
    %               .ngrid : integer specifying total number of points
    %               .r : (dim, ngrid) array containing the points
    %               .half_sidelen : float specifying the size of the regular grid.
    %         proxy_info: struct with fields
    %               .n_points_total : integer describing the total number of proxy points
    %               .dim : dimension of the problem
    %               .n_per_dim_3D : integer specifying how many points per dimension for 
    %                   discretizing the cube.
    %               .radius : float specifying the distance of proxy
    %                   surface from the center of the regular grid.
    %               .r : (dim, n_points_total) array containing the points



    dim = size(src_info.r, 1);
    nsrc = size(src_info.r, 2);
    if nargin < 5
        n_nbr = 1000;
    end
    crad = 2;

    % Get the half_sidelen and center of the points to specify the regular grid
    [Lbd, center] = bounding_box([src_info.r,targ_info.r]);
    halfside = spread_halfside([src_info.r,targ_info.r], n_nbr, crad);

    % get prototype grid for spreading
    [spread_info, proxy_info] = get_nspread_and_nproxy(kernel, dim, tol, halfside);
    
    % get total grid
    ndim = ceil(diff(Lbd, 1, 2) / spread_info.dx);
    rpad = 2;

    if dim == 2
        xx = Lbd(1,1) + (0:rpad*ndim(1)) *spread_info.dx;
        yy = Lbd(2,1) + (0:rpad*ndim(2)) *spread_info.dx;
        [X, Y] = meshgrid(xx,yy);
        rgrid = [X(:).'; Y(:).'];
    elseif dim ==3
        xx = Lbd(1,1) + (0:rpad*ndim(1)) *spread_info.dx;
        yy = Lbd(2,1) + (0:rpad*ndim(2)) *spread_info.dx;
        zz = Lbd(3,1) + (0:rpad*ndim(3)) *spread_info.dx;
        [X, Y, Z] = meshgrid(xx,yy);
        rgrid = [X(:).'; Y(:).' ; Z(:).'];
    end

    


end
