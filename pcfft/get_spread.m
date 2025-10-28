function [A_spread, sort_info] = get_spread(kern_0, kern, src_info, grid_info, proxy_info)
    % This routine returns the matrix that maps charge strengths at srcinfo.r to 
    % charge strengths on the equispaced grid.
    %
    %
    %
    % Inputs: kern_0 : 
    %           kern : 
    %           src_info : struct with fields
    %               .n_per_dim : integer describing the number of points per dimension
    %               .dim : integer specifying the dimension
    %               .ngrid : integer specifying total number of points
    %               .r : (dim, ngrid) array containing the points
    %               .half_sidelen : float specifying the size of the regular grid.
    %         proxy_info: struct with fields
    %               .n_points_total : integer describing the total number of proxy points
    %               .dim : dimension of the problem
    %               .n_per_dim_3D : integer specifying how many points per dimension for 
    %                   discretizing the cube.
    %               .radius : float specifying the distance of proxy
    %                   surface from the center of the regular grid.
    %               .r : (dim, n_points_total) array containing the points
    % Returns: A_spread : A matrix mapping charge strengths at the source points to weights on the equispaced grid. Has shape (n_equispaced_pts, n_src_pts)

    dim = proxy_info.dim;


    % Evaluate kernel src -> proxy.
    K_src_to_proxy = kern_0(src_info.r, proxy_info.r);

    % Evaluate kernel equispaced -> proxy.
    K_equi_to_proxy = kern_0(grid_info.r, proxy_info.r);

    % To find the equispaced weights c, we solve a least squares problem.
    % K_equi_to_proxy c = K_src_to_proxy * src_weights
    % So c = pinv(K_equi_to_proxy) * K_src_to_proxy * src_weights.
    A_spread = pinv(K_equi_to_proxy) * K_src_to_proxy;

    sort_info = null(1);

end