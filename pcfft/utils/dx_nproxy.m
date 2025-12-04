function [spread_info, proxy_info] = dx_nproxy(kernel, dim, tol, halfside, crad)
    
    if nargin < 5
        crad = 2;
    end
    % We intend to solve least squares problems with a 
    % rank-deficient matrix in this script.
    % warning('off', 'MATLAB:rankDeficientMatrix')


    % Initialize a collection of 100 random source points and weights inside a
    % box/cube specified by size halfside
    nsrc = 100;
    rng(0);

    % source points belong to a bin which has size bwidth = c_bwidth * halfside
    % along each side. So if c_bwidth = 2, the bin is the same size as the
    % spreading box, and if c_bwidth = 1, the bin is half the size of the
    % spreading box.

    c_bwidth = 1.0;
    bwidth = c_bwidth * halfside;

    src_pts = (rand(dim, nsrc) - 0.5) * bwidth;
    src_weights = rand(nsrc, 1) - 0.5;

    % R = max distance from origin to corner of spreading box
    R = sqrt(dim) * halfside;
    % crad = r/R is a parameter in our hands.
    % radius (little r) is crad * R;
    radius = crad * R;


    % Initialize a set of eval points that are 1.1 * radius away from center 
    % of src points
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
        nspread = 10;
    else
        nspread = 5;
    end
 
    % disp("dx_nproxy: halfside: " + num2str(halfside))

    bool_unconverged = true;

    % First, increase nspread and nproxy until we meet the error tolerance
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


        % The regular grid points span [-halfside, halfside] in each dimension
        % They are sampled so the first point is at -halfside + dx/2 and the last
        % point is at halfside - dx/2
        dx =  2 * halfside / nspread;
        % The discretization points start at -halfside + dx/2
        xx = -halfside + dx / 2 + (0:nspread - 1) * dx;
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


        % Check error between approx and exact
        errors = abs(approx_at_target(:) - target_evals(:));
        err = max(errors);

        bool_unconverged = err >= tol;

        if nspread > 100
            error("Computing proxy size didn't converge after nspread > 100")
        end

    end
    % disp("dx_nproxy: Intermediate nproxy " + int2str(nproxy));
    % disp("dx_nproxy: Intermediate nspread " + int2str(nspread));

    % Compute the kernel evals from source to proxy pts
    % We don't really need this inside this loop since proxy_pts are fixed.
    K_source_to_proxy = kernel(src_pts, proxy_pts);
    % Solve the least squares problem to find weights
    evals_at_proxy = K_source_to_proxy * src_weights;
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
        reg_pts_old = reg_pts;
        dx =  2 * halfside / nspread;
        % The discretization points start at -halfside/2 + dx/2
        xx = -halfside + dx / 2 + (0:nspread - 1) * dx;

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


        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);

        spread_weights = K_reg_to_proxy \ evals_at_proxy;

        % Eval the approximation at the eval point
        approx_at_target = kernel(reg_pts, target_pts) * spread_weights;

        % Check error between approx and exact
        errors = abs(approx_at_target(:) - target_evals(:));
        err = max(errors);
        bool_converged = err < tol;

        % disp("dx_nproxy: Down pass: Trying nspread " + int2str(nspread) + " with error " + num2str(err))

        % min nspread = 4.
        if nspread == 3
            break;
        end

    end

    % Undo the last iteration which went over error tol
    dx = dx_old;
    reg_pts = reg_pts_old;
    nspread = nspread + 1;

    % disp("dx_nproxy: After down pass, nspread " + int2str(nspread));
    nbinpts = ceil(nspread / 2);
    % disp("dx_nproxy: Intermediate nbinpts " + int2str(nbinpts));

    % If nspread // 2 == 0, we can exit early because the bin size is 
    % the same as half the spreading box size.
    if mod(nspread, 2) == 0
        % disp("dx_nproxy: Final bin width: " + num2str(nbinpts * dx))
        spread_info = struct;
        proxy_info = struct;

        spread_info.dx = dx;
        spread_info.nspread = nspread;
        spread_info.nbinpts = nbinpts;

        proxy_info.dim = dim;
        if dim == 3
            proxy_info.n_points_total = 6 * (nproxy ^ 2);
            proxy_info.n_per_dim_3D = nproxy;
        else
            proxy_info.n_points_total = nproxy;
        end
        proxy_info.r = proxy_pts;
        proxy_info.radius = radius;
        return;
    end
    
    % Otherwise, the bin may be slightly wider than half the spreading box,
    % so we will re-draw random source points in the new bin size
    bwidth = nbinpts * dx;
    % disp("dx_nproxy: Final bin width: " + num2str(bwidth))


    src_pts_newbin = (rand(dim, nsrc) - 0.5) * bwidth;
    % src_pts_newbin(:, 1:nsrc_deterministic) = src_pts_d;
    K_src_to_target_new = kernel(src_pts_newbin, target_pts);
    target_evals_new = K_src_to_target_new * src_weights(:);

    % Loop to find a final nproxy with the new bin size
    bool_unconverged = true;
    nproxy = nproxy - 1;
    while bool_unconverged
        nproxy = nproxy + 1;
        if dim == 2
            proxy_pts = get_ring_points(nproxy, radius);
        else
            proxy_pts = get_sphere_points(nproxy, radius);
        end

        % Compute the kernel evals from source to proxy pts
        K_source_to_proxy = kernel(src_pts_newbin, proxy_pts);
        % Kernel regular -> proxy pts
        K_reg_to_proxy = kernel(reg_pts, proxy_pts);
        % Solve the least squares problem to find weights
        evals_at_proxy = K_source_to_proxy * src_weights;
        spread_weights = K_reg_to_proxy \ evals_at_proxy;
        % Eval the approximation at the eval point
        approx_at_target = kernel(reg_pts, target_pts) * spread_weights;
        % Check error between approx and exact
        errors = abs(approx_at_target(:) - target_evals_new(:));
        err = max(errors);
        bool_unconverged = err >= tol;

        % disp("dx_nproxy: Trying nproxy " + int2str(nproxy) + " with error " + num2str(err))

        if nproxy > 200
            error("dx_nproxy: Computing final proxy size didn't converge after nproxy > 200")
        end

    end

    % disp("dx_nproxy: dx: "  + num2str(dx))
    % disp("dx_nproxy: Final nproxy " + int2str(nproxy));
    % disp("dx_nproxy: Final nspread " + int2str(nspread));


    spread_info = struct;
    proxy_info = struct;

    spread_info.dx = dx;
    spread_info.nspread = nspread;
    spread_info.nbinpts = nbinpts;
    

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
