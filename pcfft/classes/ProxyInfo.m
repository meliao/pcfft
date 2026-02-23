classdef ProxyInfo
    % Describe the proxy surface for spreading weights.
    %
    % Store metadata describing the proxy surface used by the spreading
    % routines to compute spreading weights. The proxy surface can consist of
    % one or more concentric shells of points.
    %
    % Attributes
    % ----------
    % dim : int
    %   problem dimension (integer).
    % n_points_total : int
    %   total number of proxy points = nproxy * nshell.
    % nproxy : int
    %   number of proxy points per shell.
    % nshell : int
    %   number of concentric shells in the proxy.
    % halfside : float
    %   TODO
    % crad : float
    %   TODO
    % tol : float
    %   relative L_inf error tolerance for generating the proxy surface.
    % radius : float
    %   proxy radius.
    % r : (dim, n_points_total) matrix
    %   array of proxy point coords.

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
    end
    methods
        function obj = ProxyInfo(dim, n_points_total, nproxy, nshell, halfside, crad, tol, radius, r)
            obj.dim = dim;
            obj.n_points_total = n_points_total;
            obj.nproxy = nproxy;
            obj.nshell = nshell;
            obj.halfside = halfside;
            obj.crad = crad;
            obj.tol = tol;
            obj.radius = radius;
            obj.r = r;
        end
    end
end
