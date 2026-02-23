classdef GridInfo
    % Describe the regular grid used for spreading.
    %
    % Store information about the regular grid used by the spreading
    % routines. The grid is divided into non-overlapping spreading bins, each
    % of which is spread to a larger spreading box independently.
    %
    % Attributes
    % ----------
    % ngrid : (dim, 1) matrix
    %   number of regular grid points per dimension.
    % Lbd : (dim, 2) matrix
    %   containing the points (xmin, ymin, zmin) and
    %   (xmax, ymax, zmax). Describes the bounding box of the union of source
    %   and target points.
    % dx : float
    %   grid spacing.
    % nspread : int
    %   number of spreading points per dimension required for the requested
    %   accuracy.
    % nbinpts : int
    %   width in dx's of the spreading bin.
    % rpad : int
    %   TODO
    % r : (dim, npts) matrix
    %   padded regular grid points.
    % dim : int
    %   problem dimension.
    % nbin : (dim, 1) matrix
    %   number of spreading bins per dimension.
    % offset : int
    %   Number of dx's used to pad the grid.
    % rmax : (dim, 1) matrix
    %   maximum coordinate value of the padded grid.
    % rmin : (dim, 1) matrix
    %   minimum coordinate value of the padded grid.
    % n_nbr : int
    %   average number of near-field neighbours.

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
        rmax
        rmin
        n_nbr
    end
    methods
        function obj = GridInfo(Lbd, dx, nspread, nbinpts, dim, n_nbr)


    
            bin_sidelen = dx * nbinpts;

            % Number of spreading bins in each dimension.
            n_bin = ceil(diff(Lbd, 1, 2) / bin_sidelen);

            % Number of points padding each side
            pad = ceil((nspread - nbinpts) / 2);
            % Width below the bottom corner of Lbd to start the regular grid points
            offset = pad * dx - dx / 2;

            ngrid = n_bin * nbinpts + pad * 2;

            if dim == 2
                % Create a regular grid with spacing dx starting at the xmin, ymin point
                % specified by Lbd. 
                xx = Lbd(1, 1) - offset + (0: ngrid(1) - 1) * dx;
                yy = Lbd(2, 1) - offset + (0: ngrid(2) - 1) * dx;
                [X, Y] = meshgrid(xx, yy);
                rgrid = [X(:).'; Y(:).'];
            elseif dim == 3
                xx = Lbd(1, 1) - offset + (0: ngrid(1) - 1) * dx;
                yy = Lbd(2, 1) - offset + (0: ngrid(2) - 1) * dx;
                zz = Lbd(3, 1) - offset + (0: ngrid(3) - 1) * dx;
                [X, Y, Z] = meshgrid(xx, yy, zz);
                X = permute(X,[3,1,2]);
                Y = permute(Y,[3,1,2]);
                Z = permute(Z,[3,1,2]);
                rgrid = [X(:).'; Y(:).'; Z(:).'];
            end

            rmax = max(rgrid, [], 2);
            rmin = min(rgrid, [], 2);

            obj.ngrid = ngrid;
            obj.Lbd = Lbd;
            obj.dx = dx;
            obj.nspread = nspread;
            obj.nbinpts = nbinpts;
            obj.rpad = pad;
            obj.r = rgrid;
            obj.dim = dim;
            obj.nbin = n_bin;
            obj.offset = offset;
            obj.rmax = rmax;
            obj.rmin = rmin;
            obj.n_nbr = n_nbr;
        end
    end
end
