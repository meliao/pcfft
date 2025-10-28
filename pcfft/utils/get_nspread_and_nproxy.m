function [grid_info, proxy_info] = get_nspread_and_nproxy(kernel, dim, tol, halfside)
    

    % Initialize a collection of 100 random source points and weights inside a box/cube
    % specified by size halfside
    nsrc = 100;
    rng(0);
    src_pts = (rand(dim, nsrc) - 0.5) * halfside;
    src_weights = rand(nsrc, 1) - 0.5;

    % R = max distance from origin of one of our source points
    R = sqrt(dim) * halfside;
    % c = r/R is a parameter in our hands.
    c = 2;
    % radius (little r) is c * R;
    radius = c * R;


    % Initialize a set of eval points that are 1.1 * radius away from center of src points
    ntarget = 100;
    if dim == 2
        target_pts = get_ring_points(ntarget, 1.1 * radius);
    else
        % ntarget_eff = floor((ntarget / 6)^(1/3));
        target_pts = get_sphere_points(ntarget, 1.1 * radius);
    end
    K_src_to_target = kernel(src_pts, target_pts);
    target_evals = K_src_to_target * src_weights(:);


    % Start with small numbers for nspread and nproxy
    if dim == 2
        nspread = 50;
    else
        nspread = 10;
    end
 

    bool_unconverged = true;

    while bool_unconverged

        if dim == 2
            nspread = nspread + 2;
            nproxy = nspread;
        else
            nspread = nspread + 1;
            nproxy = 2 * nspread^2;
        end

        disp("Evaluating with nspread  = " + int2str(nspread))
        disp("Evaluating with nproxy = " + int2str(nproxy))

        % The "proxy" points are the ones at which we try to match the
        % outputs.

        if dim == 2
            proxy_pts = get_ring_points(nproxy, radius);
        else
            proxy_pts = get_sphere_points(nproxy, radius);
        end

        % Discretize the regular grid using nspread points
        reg_pts = get_regular_grid(nspread, halfside, dim);

        % Compute the kernel evals from source to proxy pts
        K_source_to_proxy = kernel(src_pts, proxy_pts);
        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);

        % Solve the least squares problem to find weights
        evals_at_proxy = K_source_to_proxy * src_weights;

        % warning('off', 'MATLAB:rankDeficientMatrix')
        spread_weights = K_reg_to_proxy \ evals_at_proxy;

        % Eval the approximation at the eval point
        approx_at_target = kernel(reg_pts, target_pts) * spread_weights;
        % disp("approx_at_target shape")
        % disp(size(approx_at_target));
        % disp("target_evals");
        % disp(size(target_evals));

        % Check error between approx and exact
        errors = abs(approx_at_target(:) - target_evals(:));
        err = max(errors);
        disp("Error: ");
        disp(err);

        bool_unconverged = err >= tol;

        if nspread > 500
            error("Computing proxy size didn't converge after nspread > 500")
        end

    end
    disp("Intermediate nproxy " + int2str(nproxy));
    disp("Intermediate nspread " + int2str(nspread));


    % Keep nproxy fixed and reduced nspread while keeping accuracy
    bool_converged = true;
    while bool_converged

        nspread = nspread - 1;
        % disp("Evaluating with nspread  = " + int2str(nspread))
        % disp("Evaluating with nproxy = " + int2str(nproxy))

        % The "proxy" points are the ones at which we try to match the
        % outputs.

        if dim == 2
            proxy_pts = get_ring_points(nproxy, radius);
        else
            proxy_pts = get_sphere_points(nproxy, radius);
        end

        % Discretize the regular grid using nspread points
        reg_pts = get_regular_grid(nspread, halfside, dim);

        % Compute the kernel evals from source to proxy pts
        K_source_to_proxy = kernel(src_pts, proxy_pts);
        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);

        % Solve the least squares problem to find weights
        evals_at_proxy = K_source_to_proxy * src_weights;

        % warning('off', 'MATLAB:rankDeficientMatrix')
        spread_weights = K_reg_to_proxy \ evals_at_proxy;

        % Eval the approximation at the eval point
        approx_at_target = kernel(reg_pts, target_pts) * spread_weights;
        % disp("approx_at_target shape")
        % disp(size(approx_at_target));
        % disp("target_evals");
        % disp(size(target_evals));

        % Check error between approx and exact
        errors = abs(approx_at_target(:) - target_evals(:));
        err = max(errors);
        % disp("Error: ");
        % disp(err);

        bool_converged = err < tol;


    end
    nspread = nspread + 1;

    disp("Final nproxy " + int2str(nproxy));
    disp("Final nspread " + int2str(nspread));


    grid_info = struct;
    grid_info.nspread = nspread;
    grid_info.dim = dim;
    grid_info.r = reg_pts;
    grid_info.ngrid = size(reg_pts, 2);
    grid_info.halfside = halfside;
    

    proxy_info = struct;
    proxy_info.dim = dim;
    if dim == 3
        proxy_info.n_points_total = 6 * (nproxy ^ 2);
        proxy_info.n_per_dim_3D = nproxy;
    else
        proxy_info.n_points_total = nproxy;
    end
    proxy_info.r = proxy_pts;
    proxy_info.radius = radius;

end
