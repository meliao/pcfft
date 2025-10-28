function [grid_info, proxy_info] = get_grid(kernel, src_info, targ_info, tol, fill)
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
    %           fill : (optional) int specifying the max number of source points that are allowed 
    %               to be within 3 * dx of any source point. Defaults to floor(n_src_pts / 20)
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
        fill = max(1, floor(nsrc / 20));
    end

    % Get the half_sidelen and center of the points to specify the regular grid
    [half_sidelen, center] = bounding_box(src_info.r);

    % TODO: figure out 

    if dim == 2
        eval_pt = get_ring_points(100, 1.1 * targ_info.radius);
    else
        % 6 * 16 = 96 eval points
        eval_pt = get_cube_points(4, 1.1 * targ_info.radius);
    end
    % Compute the kernel evaluation at the eval point
    K_src_to_eval = kernel(src_info.r, eval_pt);
    k_evals = K_src_to_eval * src_info.weights;



    % Start the equispaced discretization large enough so that the fill
    % parameter is satisfied.
    n_reg_pts = get_reg_pts_satisfying_fill(src_info.r, half_sidelen, fill);

    % Starting discretization
    n_proxy_pts = 8;
    n_reg_pts = 8;

    if dim == 3
        % In 3D, count the # of proxy points via number of points
        % per dimension on each face.
        % There are 6 * (n_proxy_pts)^2 points on the cube.
        n_proxy_pts = 3;
        n_reg_pts = 3;
    end



    bool_unconverged = true;

    while bool_unconverged

        if dim == 2
            n_reg_pts = n_reg_pts + 2;
            n_proxy_pts = n_proxy_pts + 2;
        else
            n_reg_pts = n_reg_pts + 1;
            n_proxy_pts = n_proxy_pts + 1;
        end

        disp("Evaluating with n_reg_pts = " + int2str(n_reg_pts))
        disp("Evaluating with n_proxy_pts = " + int2str(n_proxy_pts))

        % The "proxy" points are the ones at which we try to match the
        % outputs.

        if dim == 2
            % Discretize the ring using n_ring_pts
            proxy_pts = get_ring_points(n_proxy_pts, targ_info.radius);
        else
            proxy_pts = get_cube_points(n_proxy_pts, targ_info.radius);
        end

        % Discretize the regular grid using n_reg_pts
        reg_pts = get_regular_grid(n_reg_pts, half_sidelen, dim);

        % Compute the kernel evals from source to proxy pts
        K_source_to_proxy = kernel(src_info.r, proxy_pts);
        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);

        % Solve the least squares problem to find weights
        evals_at_proxy = K_source_to_proxy * src_info.weights;

        % warning('off', 'MATLAB:rankDeficientMatrix')
        disp("Least squares problem size:")
        disp(size(K_reg_to_proxy))
        reg_weights = K_reg_to_proxy \ evals_at_proxy;

        % Eval the approximation at the eval point
        approx_at_eval_pt = kernel(reg_pts, eval_pt) * reg_weights;

        % Check error between approx and exact
        errors = abs(approx_at_eval_pt(:) - k_evals(:));
        err = max(errors);

        bool_unconverged = err >= tol;

        if n_reg_pts > 100
            error("Computing proxy size didn't converge after 100 points")
        end

    end

    grid_info = struct;
    grid_info.n_per_dim = n_reg_pts;
    grid_info.dim = dim;
    grid_info.r = reg_pts;
    grid_info.ngrid = size(reg_pts, 2);
    grid_info.half_sidelen = half_sidelen;
    
    proxy_info = struct;
    proxy_info.dim = dim;
    if dim == 3
        proxy_info.n_points_total = 6 * (n_proxy_pts ^ 2);
        proxy_info.n_per_dim_3D = n_proxy_pts;
    else
        proxy_info.n_points_total = n_proxy_pts;
    end
    proxy_info.r = proxy_pts;
    proxy_info.radius = targ_info.radius;

end
