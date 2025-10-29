function [A_spread, sort_info] = get_spread(kern_0, kern, src_info, grid_info, proxy_info, nbin)
    % This routine returns the matrix that maps charge strengths at srcinfo.r to 
    % charge strengths on the equispaced grid.
    % We want to 
    if nargin < 6
        nbin = 1;
    end
    dim = proxy_info.dim;


    % First, sort the points into bins
    if dim == 2
        [r_sorted, bin_idxes, id_start] = bin_pts_2d(src_info.r, grid_info.dx, grid_info.ngrid, grid_info.Lbd, nbin)
    else
        [r_sorted, bin_idxes, id_start] = bin_pts_3d(src_info.r, grid_info.dx, grid_info.ngrid, grid_info.Lbd, nbin)
    end

    % Now, loop through the bins and 


end