function [n_reg_pts, n_proxy_pts] = find_proxy_size(kernel, ...
        half_sidelen, dim, radius, tol)
    % Computes the number of regular gridpoints on a box or cube with size
    % <half_sidelen> needed to approximate an evaluation of the kernel at a
    % point at least <radius> away.
    %
    %
    % kernel is a function handle accepting (source_pts, target_pts)
    % half_sidelen float describing the size of the square or cube
    % dim integer specifying the dimension of the domain
    % radius is a float
    % tol is a float
    %
    % Returns n_reg_pts, the number of regular grid points in one dimension
    %         n_distant_pts, the number of points on a distant ring (2D) or
    %         cube (3D).

    rng(0)
    n_src = 100;
    src_pts = (rand(dim, n_src) - 0.5) * half_sidelen;
    src_weights = rand(n_src);

    if dim == 2
        eval_pt = get_ring_points(100, 1.1 * radius);
    else
        eval_pt = get_cube_points(100, 1.1 * radius);
    end

    % Starting discretization
    n_proxy_pts = 8;
    n_reg_pts = 8;
    disp(dim);

    % if dim == 3
    %     % In 3D, count the # of proxy points via number of points
    %     % per dimension on each face.
    %     % There are 6 * (n_proxy_pts)^2 points on the cube.
    %     n_proxy_pts = 3;
    % end

    % Compute the kernel evaluation at the eval point
    K_src_to_eval = kernel(src_pts, eval_pt);
    k_evals = K_src_to_eval * src_weights;

    bool_unconverged = true;

    while bool_unconverged

        n_reg_pts = n_reg_pts + 2;
        n_proxy_pts = n_proxy_pts + 2;

        % The "proxy" points are the ones at which we try to match the
        % outputs.

        if dim == 2
            % Discretize the ring using n_ring_pts
            proxy_pts = get_ring_points(n_proxy_pts, radius);
        else
            proxy_pts = get_cube_points(n_proxy_pts, radius);
        end

        % Discretize the regular grid using n_reg_pts
        reg_pts = get_regular_grid(n_reg_pts, half_sidelen, dim);

        % Compute the kernel evals from source to proxy pts
        K_source_to_proxy = kernel(src_pts, proxy_pts);
        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);

        % Solve the least squares problem to find weights
        evals_at_proxy = K_source_to_proxy * src_weights;

        % warning('off', 'MATLAB:rankDeficientMatrix')
        reg_weights = K_reg_to_proxy \ evals_at_proxy;

        % Eval the approximation at the eval point
        approx_at_eval_pt = kernel(reg_pts, eval_pt) * reg_weights;

        % Check error between approx and exact
        errors = abs(approx_at_eval_pt(:) - k_evals(:));
        disp("Shape of errors:")
        disp(size(errors))
        err = max(errors);

        bool_unconverged = err >= tol;

        if n_reg_pts > 100
            error("Computing proxy size didn't converge after 100 points")
        end

    end

    if dim == 3
        n_proxy_pts = 6 * (n_proxy_pts ^ 2);
    end

end
