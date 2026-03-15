function [grid_info, proxy_info] = get_grid(kernel, src_info, targ_info, ...
        tol, n_nbr, opts)
    % Compute the equispaced spreading grid used for the precorrected FFT.
    %
    % Parameters
    % ----------
    % kernel : kernel
    %   Free-space kernel.
    % src_info : point_info
    %   struct giving information about the sources.
    % targ_info : point_info
    %   struct describing the target points.
    % tol : float
    %   float specifying absolute error tolerance. Error is evaluated at a
    %   surface 1.1 * radius of the innermost proxy surface.
    % n_nbr : int, optional
    %   int specifying the average number of interactions that must be done
    %   directly. Defaults to 1000.
    % opts : struct, optional
    %   options to manipulate the choice of proxy points. Available
    %   options:
    %
    %   - opts.multi_shells - Whether to default to shell-based proxies. 
    %           Defaults to false. Accelerates the precomputation for
    %           kernels that are known to not satisfy Green's identity
    %   - opts.proxy_der - Number of radial derivatives to use in the proxy 
    %           Can be a number between 0 and 2. Defaults to 0. This option
    %           can avoid invoking shells for kernels that are derived from
    %           high order PDEs. (see wrap_kern_der)
    %   - opts.halfside - Manually set box size (and implicitly n_nbr).
    %           Only recommended for expert users. (See pcff_fmm3dbie_demo.m)
    %           Useful for plotting BIE solutions without overly small boxes.
    %
    %
    % Returns
    % -------
    % grid_info : GridInfo
    %   GridInfo object specifying the regular grid used for spreading.
    % proxy_info : ProxyInfo
    %   ProxyInfo object specifying the proxy surface(s) used in the proxy point method.


    dim = size(src_info.r(:,:), 1);
    if nargin < 5 || isempty(n_nbr)
        n_nbr = 1000;
    end
    if nargin < 6
        opts = false;
    end
    multi_shells = false;
    if isfield(opts,'multi_shells')
        multi_shells = opts.multi_shells;
    end
    proxy_der = 0;
    if isfield(opts,'proxy_der')
        proxy_der = opts.proxy_der;
    end

    if ~isa(kernel,'function_handle')
        try
            kernel = kernel.eval;
        catch
            error('kernel is not a function and does not have an eval property')
        end
    end

    crad = 2;

    % Get the half_sidelen and center of the points to specify the regular grid
    [Lbd, ~] = bounding_box([src_info.r(:,:), targ_info.r(:,:)]);
    if isfield(opts,'halfside')
        halfside = opts.halfside;
    else
        halfside = spread_halfside([src_info.r(:,:), targ_info.r(:,:)], n_nbr, crad);
    end
    % get prototype grid for spreading
    [dx, nspread, nbinpts, proxy_info] = dx_nproxy(kernel, dim, tol, halfside, crad, multi_shells, proxy_der);


    grid_info = GridInfo(Lbd, dx, nspread, nbinpts, dim, n_nbr);
    
end
