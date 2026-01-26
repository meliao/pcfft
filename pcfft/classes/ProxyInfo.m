classdef ProxyInfo
    %PROXYINFO Describes the proxy surface for spreading weights.
    %
    %   PROXYINFO stores metadata describing the proxy surface used by the
    %   spreading routines to compute spreading weights. The proxy surface can
    %   consist of one or more concentric shells of points.
    %
    %   Properties
    %       dim              - problem dimension (integer)
    %       n_points_total   - total number of proxy points = nproxy * nshell
    %       nproxy           - number of proxy points per shell
    %       nshell           - number of concentric shells in the proxy
    %       halfside         - TODO
    %       crad             - TODO
    %       tol              - relative L_inf error tolerance for generating 
    %                          the proxy surface
    %       radius           - proxy radius
    %       r                - (dim, n_points_total) array of proxy point coords

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