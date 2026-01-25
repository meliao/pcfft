classdef GridInfo
    %GRIDINFO Describes the regular grid used for spreading.
    %
    %   GRIDINFO stores information about the regular grid used by the
    %   spreading routines. The grid is divided into non-overlapping spreading
    %   bins, each of which is spread to a larger spreading box independently.
    %
    %   Properties
    %       ngrid    - (dim,1) number of regular grid points per dimension
    %       Lbd      - (dim,2) array containing the points 
    %                       (xmin, ymin, zmin) and (xmax, ymax, zmax). Describes
    %                       the bounding box of the union of source and target
    %                       points. 
    %       dx       - grid spacing 
    %       nspread  - number of spreading points per dimension required
    %                  for the requested accuracy
    %       nbinpts  - width in dx's of the spreading bin
    %       rpad     - TODO
    %       r        - (dim, npts) array of padded regular grid points
    %       dim      - problem dimension (integer)
    %       nbin     - (dim,1) number of spreading bins per dimension
    %       offset   - (integer) offset used for indexing the padded grid

    properties
        ngrid
        Lbd
        dx
        nspread
        nbinpts
        rpad
        r
        dim
        nbin
        offset
    end
    methods
        function obj = GridInfo(ngrid, Lbd, dx, nspread, nbinpts, rpad, r, dim, nbin, offset)
            obj.ngrid = ngrid;
            obj.Lbd = Lbd;
            obj.dx = dx;
            obj.nspread = nspread;
            obj.nbinpts = nbinpts;
            obj.rpad = rpad;
            obj.r = r;
            obj.dim = dim;
            obj.nbin = nbin;
            obj.offset = offset;
        end
    end
end