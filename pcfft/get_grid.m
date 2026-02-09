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
    % Returns grid_info: GridInfo object
    %         proxy_info: ProxyInfo object

    dim = size(src_info.r(:,:), 1);
    % nsrc = size(src_info.r, 2);
    if nargin < 5
        n_nbr = 1000;
    end

    crad = 2;

    % Get the half_sidelen and center of the points to specify the regular grid
    [Lbd, ~] = bounding_box([src_info.r(:,:), targ_info.r(:,:)]);
    halfside = spread_halfside([src_info.r(:,:), targ_info.r(:,:)], n_nbr, crad);

    % get prototype grid for spreading
    [dx, nspread, nbinpts, proxy_info] = dx_nproxy(kernel, dim, tol, halfside, crad);


    grid_info = GridInfo(Lbd, dx, nspread, nbinpts, dim, n_nbr);
    
    bin_sidelen = dx * nbinpts;

    % Number of spreading bins in each dimension.
    n_bin = ceil(diff(Lbd, 1, 2) / bin_sidelen);

    % Number of points padding each side
    pad = ceil((nspread - nbinpts) / 2);
    % Width below the bottom corner of Lbd to start the regular grid points
    offset = pad * dx - dx / 2;

    ngrid = n_bin * nbinpts + pad * 2;

    disp("get_grid: ngrid:")
    disp(ngrid)

    if dim == 2
        % Create a regular grid with spacing dx starting at the xmin, ymin point
        % specified by Lbd. 
        xx = Lbd(1, 1) - offset + (0: ngrid(1) - 1) * dx;
        yy = Lbd(2, 1) - offset + (0: ngrid(2) - 1) * dx;
        [X, Y] = meshgrid(xx, yy);
        rgrid = [X(:).'; Y(:).'];
    elseif dim == 3
        xx = Lbd(1, 1) - offset + (0: ngrid(1) - 1) * dx;
        yy = Lbd(2, 1) - offset + (0: ngrid(2) - 1) * dx;
        zz = Lbd(3, 1) - offset + (0: ngrid(3) - 1) * dx;
        [X, Y, Z] = meshgrid(xx, yy, zz);
        X = permute(X,[3,1,2]);
        Y = permute(Y,[3,1,2]);
        Z = permute(Z,[3,1,2]);
        rgrid = [X(:).'; Y(:).'; Z(:).'];
    end

    % xx_dx = xx(2) - xx(1);
    % disp("get_grid: dx of xx: " + num2str(xx_dx))

    % yy_dy = yy(2) - yy(1);
    % disp("get_grid: dy of yy: " + num2str(yy_dy))
end
