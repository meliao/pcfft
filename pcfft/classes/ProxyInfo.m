classdef ProxyInfo
    % Describes the proxy surface for spreading weights.
    %
    % Stores metadata describing the proxy surface used by the spreading
    % routines to compute spreading weights. The proxy surface can consist of
    % one or more concentric shells of points.
    %
    % Attributes
    % ----------
    % dim : int
    %   problem dimension (integer).
    % n_points_total : int
    %   total number of proxy points = ``nproxy`` * ``nshell``.
    % nproxy : int
    %   number of proxy points per shell.
    % nshell : int
    %   number of concentric shells in the proxy.
    % halfside : float
    %   Specifies the halfside of the spreading box.
    % crad : float
    %   Parameter for determining the radius of the innermost proxy shell.
    % tol : float
    %   relative L_inf error tolerance for generating the proxy surface.
    % radius : float
    %   Radius of the innermost proxy shell. Used for determining near-field neighbors. = ``sqrt(dim)`` * ``halfside`` * ``crad``.
    % r : array [dim, n_points_total]
    %   array of proxy point coords for a spreading box  centered at the origin.
    % proxy_der : int
    %   number of derivatives used in proxy compression. See wrap_kern_der.

    properties
        dim
        n_points_total
        nproxy
        nshell
        halfside
        crad
        tol
        radius
        r
        proxy_der
    end
    methods
        function obj = ProxyInfo(dim, n_points_total, nproxy, nshell, halfside, crad, tol, radius, r, proxy_der)
            obj.dim = dim;
            obj.n_points_total = n_points_total;
            obj.nproxy = nproxy;
            obj.nshell = nshell;
            obj.halfside = halfside;
            obj.crad = crad;
            obj.tol = tol;
            obj.radius = radius;
            obj.r = r;
            obj.proxy_der = proxy_der;
        end
    end
end
