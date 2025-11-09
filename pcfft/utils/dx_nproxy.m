function [spread_info, proxy_info] = dx_nproxy(kernel, dim, tol, halfside, crad)
    
    if nargin < 5
        crad = 2;
    end
    % We intend to solve least squares problems with a 
    % rank-deficient matrix in this script.
    warning('off', 'MATLAB:rankDeficientMatrix')


    % Initialize a collection of 100 random source points and weights inside a box/cube
    % specified by size halfside
    nsrc = 100;
    rng(0);
    src_pts = (rand(dim, nsrc) - 0.5) * halfside;
    src_weights = rand(nsrc, 1) - 0.5;

    % disp("dx_nproxy: src_pts max:")
    % disp(max(src_pts))

    % R = max distance from origin of one of our source points
    R = sqrt(dim) * halfside;
    % crad = r/R is a parameter in our hands.
    % radius (little r) is crad * R;
    radius = crad * R;


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

        % disp("Evaluating with nspread  = " + int2str(nspread))
        % disp("Evaluating with nproxy = " + int2str(nproxy))

        % The "proxy" points are the ones at which we try to match the
        % outputs.

        if dim == 2
            proxy_pts = get_ring_points(nproxy, radius);
        else
            proxy_pts = get_sphere_points(nproxy, radius);
        end

        % Discretize the box [-halfside / 2, halfside / 2]^dim
        % using nspread points per dimension.
        reg_pts = get_regular_grid(nspread, halfside / 2, dim);

        % Compute the kernel evals from source to proxy pts
        K_source_to_proxy = kernel(src_pts, proxy_pts);
        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);

        % Solve the least squares problem to find weights
        evals_at_proxy = K_source_to_proxy * src_weights;

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
    disp("dx_nproxy: Intermediate nproxy " + int2str(nproxy));
    disp("dx_nproxy: Intermediate nspread " + int2str(nspread));


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
        reg_pts_old = reg_pts;
        reg_pts = get_regular_grid(nspread, halfside / 2, dim);

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


    % First, compute dx
    dx = reg_pts_old(2, 2) - reg_pts_old(2, 1);
    disp("dx_nproxy: dx: "  + num2str(dx))

    % Now, nspread = number of grid points with spacing dx needed to cover 
    % the interval [- halfside/2, halfside/2]
    nspread = ceil(halfside / dx);

    % The regular points now span a box with size 
    % dx * (nspread + 1).
    % This can be slightly larger than the original halfside
    new_halfside = dx * (nspread + 1);
    src_pts_2 = (rand(dim, nsrc) - 0.5) * new_halfside;
    % disp("dx_nproxy: src_pts_2 max: ")
    % disp(max(src_pts_2))
    K_src_to_target = kernel(src_pts_2, target_pts);
    target_evals = K_src_to_target * src_weights(:);

    % We want to compute a new set of (dx, nproxy, nspread)
    % params for this spreading box width.

    bool_unconverged = true;

    while bool_unconverged

        if dim == 2
            nspread = nspread + 2;
            nproxy = nproxy  + 2;
            proxy_pts = get_ring_points(nproxy, radius);

        else
            nspread = nspread + 1;
            nproxy = 2 * nspread^2;
            proxy_pts = get_sphere_points(nproxy, radius);

        end


        % Discretize the box [-new_halfside / 2, new_halfside / 2]^dim
        % using nspread points per dimension.
        reg_pts_old = reg_pts;
        reg_pts = get_regular_grid(nspread, new_halfside /2, dim);

        % Compute the kernel evals from source to proxy pts
        K_source_to_proxy = kernel(src_pts_2, proxy_pts);
        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);

        % Solve the least squares problem to find weights
        evals_at_proxy = K_source_to_proxy * src_weights;

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

        if nspread > 50
            error("Computing proxy size didn't converge after nspread > 500")
        end

    end


    disp("dx_nproxy: Final nproxy " + int2str(nproxy));
    disp("dx_nproxy: Final nspread " + int2str(nspread));

    % First, compute dx
    dx = reg_pts_old(2, 2) - reg_pts_old(2, 1);
    disp("dx_nproxy: dx: "  + num2str(dx))
    nspread = ceil(new_halfside / dx);


    spread_info = struct;
    proxy_info = struct;

    spread_info.dx = dx;
    spread_info.nspread = nspread;
    % spread_info.dim = dim;
    % spread_info.rloc = reg_pts;
    % spread_info.ngrid = size(reg_pts, 2);
    % spread_info.halfside = halfside;
    

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
