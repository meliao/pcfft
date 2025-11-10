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
    nsrc_deterministic = 2^dim ;
    % Put source points at the corner of the spreading grid.
    if dim == 2
        src_pts_d = [-halfside/2 -halfside/2
                    -halfside/2 halfside/2
                    halfside/2 -halfside/2
                    halfside/2 halfside/2].';
    else
        src_pts_d = [-halfside/2 -halfside/2 -halfside/2
                    -halfside/2 halfside/2 -halfside/2
                    halfside/2 -halfside/2 -halfside/2
                    halfside/2 halfside/2  -halfside/2
                    -halfside/2 -halfside/2 halfside/2
                    -halfside/2 halfside/2 halfside/2
                    halfside/2 -halfside/2 halfside/2
                    halfside/2 halfside/2  halfside/2].';
    end

    % src_pts = [-halfside/2 halfside/2] * ones(2, dim);
    % assert(all(size(src_pts) == [2 dim]));
    src_pts = (rand(dim, nsrc) - 0.5) * halfside;
    src_pts(:, 1:nsrc_deterministic) = src_pts_d;
    src_weights = rand(nsrc, 1) - 0.5;
    disp("src_pts shape: ")
    disp(size(src_pts))

    % disp("dx_nproxy: src_pts max:")
    % disp(max(src_pts))

    % R = max distance from origin of one of our source points
    R = sqrt(dim) * halfside;
    % crad = r/R is a parameter in our hands.
    % radius (little r) is crad * R;
    radius = crad * R;


    % Initialize a set of eval points that are  radius away from center of src points
    ntarget = 100;
    if dim == 2
        target_pts = get_ring_points(ntarget, radius);
    else
        % ntarget_eff = floor((ntarget / 6)^(1/3));
        target_pts = get_sphere_points(ntarget, radius);
    end
    K_src_to_target = kernel(src_pts, target_pts);
    target_evals = K_src_to_target * src_weights(:);


    % Start with small numbers for nspread and nproxy
    if dim == 2
        nspread = 50;
    else
        nspread = 10;
    end
 
    disp("dx_nproxy: halfside: " + num2str(halfside))

    bool_unconverged = true;

    while bool_unconverged

        if dim == 2
            nspread = nspread + 2;
            nproxy = nspread;
            proxy_pts = get_ring_points(nproxy, radius);

        else
            nspread = nspread + 1;
            nproxy = 2 * nspread^2;
            proxy_pts = get_sphere_points(nproxy, radius);

        end


        dx =  halfside / nspread;
        % disp("dx_nproxy: dx: " + num2str(dx))

        % The discretization points start at -halfisde/2 + dx/2
        xx = -halfside/2 + dx / 2 + (0:nspread - 1) * dx;
        % disp("dx_nproxy: xx")
        % disp(xx)
        yy = xx;

        if dim == 2
            [X, Y] = meshgrid(xx, yy);
            reg_pts = [X(:).'; Y(:).'];
        else
            zz = xx;
            [X, Y, Z] = meshgrid(xx, yy, zz);
            X = permute(X,[3,1,2]);
            Y = permute(Y,[3,1,2]);
            Z = permute(Z,[3,1,2]);
            reg_pts = [X(:).'; Y(:).'; Z(:).'];
        end


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
        % disp("Error: ");
        % disp(err);

        bool_unconverged = err >= tol;

        if nspread > 100
            error("Computing proxy size didn't converge after nspread > 100")
        end

    end
    disp("dx_nproxy: Intermediate nproxy " + int2str(nproxy));
    disp("dx_nproxy: Intermediate nspread " + int2str(nspread));


    % Keep nproxy fixed and reduced nspread while keeping accuracy
    bool_converged = true;
    while bool_converged

        nspread = nspread - 1;


        if dim == 2
            proxy_pts = get_ring_points(nproxy, radius);
        else
            proxy_pts = get_sphere_points(nproxy, radius);
        end


        dx_old = dx;
        dx =  halfside / nspread;
        % disp("dx_nproxy: dx: " + num2str(dx))

        % The discretization points start at -halfisde/2 + dx/2
        xx = -halfside/2 + dx / 2 + (0:nspread - 1) * dx;
        % disp("dx_nproxy: xx")
        % disp(xx)
        yy = xx;

        if dim == 2
            [X, Y] = meshgrid(xx, yy);
            reg_pts = [X(:).'; Y(:).'];
        else
            zz = xx;
            [X, Y, Z] = meshgrid(xx, yy, zz);
            X = permute(X,[3,1,2]);
            Y = permute(Y,[3,1,2]);
            Z = permute(Z,[3,1,2]);
            reg_pts = [X(:).'; Y(:).'; Z(:).'];
        end

        % Discretize the regular grid using nspread points


        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);


        % In this while loop, proxy points are fixed so we don't need to 
        % recompute K_source_to_proxy
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

        bool_converged = err < tol;

    end

    % Undo the last iteration which went over error tol
    dx = dx_old;
    nspread = nspread + 1;

    disp("dx_nproxy: dx: "  + num2str(dx))
    disp("dx_nproxy: Final nproxy " + int2str(nproxy));
    disp("dx_nproxy: Final nspread " + int2str(nspread));

    % First, compute dx
    disp("dx_nproxy: dx: "  + num2str(dx))


    spread_info = struct;
    proxy_info = struct;

    spread_info.dx = dx;
    spread_info.nspread = nspread;
    

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
