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
    %   - opts.multi_shells
    %           Whether to allow more than one proxy shell. Defaults to true.
    %   - opts.proxy_der
    %           Number of radial derivatives to use in the proxy. Can be a
    %           number between 0 and 2. Defaults to 0. This option can
    %           reduce the total number of proxy points for higher order
    %           PDE problems. (see wrap_kern_der)
    %   - opts.halfside
    %           Manually set box size (and implicitly n_nbr).
    %           Only recommended for expert users. (See
    %           pcff_fmm3dbie_demo.m) Useful for plotting BIE solutions
    %           where halfside should be set by the boundary.
    %   - opts.nproxy_max
    %           Maximum number of proxy points per shell (default 6000)
    %   - opts.nshell_max
    %           Maximum number of shells (default 20)
    %   - opts.nspread_max
    %           Maximum grid points across the box (default 100 in 2D, 30 in 3D)
    %   - opts.relax_nshell
    %           Trade one more nspread for a smaller nshell (default false)
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
    if nargin < 6 || ~isstruct(opts)
        opts = struct();
    end
    multi_shells = true;
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
    [dx, nspread, nbinpts, proxy_info] = dx_nproxy(kernel, dim, tol, halfside, ...
        crad, multi_shells, proxy_der, opts);


    grid_info = GridInfo(Lbd, dx, nspread, nbinpts, dim, n_nbr, proxy_info.radius);
    
end
