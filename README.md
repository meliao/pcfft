# pcfft
This repository provides routines for quickly computing N-Body calculations in 2D and 3D and is designed to work with generic translation invariant kernels



### User-callable precomputation routines
[grid_info, pxyinfo] = get_grid(kern, srcinfo, targinfo, eps, fill)
This routine determines the size of the equispaced grid, the number of proxy points used for spreading, and the spreading parameters

[A_spread, sortinfo] = get_spread(kern_0, kern, srcinfo, grid_info, pxyinfo)
This routine returns the matrix that maps charge strengths at srcinfo.r to charge strengths on the equispaced grid, it also returns some point binning info used in get_addsub

A_add_sub = get_addsub(kern_0, kern_s, kern_t, kern_st, srcinfo, targinfo, grid_info, pxyinfo, sortinfo_s, sortinfo_t)







