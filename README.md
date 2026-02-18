# pcfft
This repository provides routines for quickly computing N-Body calculations in 2D and 3D and is designed to work with generic translation invariant kernels

$$u_i = \sum_{j\neq i} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j$$

The routines can also be used to compute sums involving source and target derivatives, such as

$$w_i = \sum_{j\neq i} \partial_{\mathbf{n}_i}\partial_{\mathbf{n}_j} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j$$

### User-callable precomputation routines
[grid_info, pxy_info] = get_grid(kern_0, srcinfo, targinfo, eps, fill)
This routine determines the size of the equispaced grid, the number of proxy points used for spreading, and the spreading parameters

[A_spread, K_src_to_reg, sort_info] = get_spread(kern_0, kern_der, srcinfo, grid_info, pxy_info)
This routine returns the matrix that maps charge strengths at srcinfo.r to charge strengths on the equispaced grid, it also returns some point binning info used in get_addsub

A_add_sub = get_addsub(kern_0, kern_st, srcs, targs, grid_info, proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t)
This routine returns the matrix that fixes the interactions that are done incorrectly by the fast apple

kern_hat = get_kernhat(kern_0, grid_info)
evaluate FFT of kern_0 on grid


### Sample use

See the `demos/` directory for a variety of examples, including demonstrating pairing this package with some popular repositories.


### Depenencies

 * [FLAM](https://github.com/klho/FLAM)


### wishlist:
* vector valued kernels


### proposed parameter plan:

* Determine spreading box half side length
* Sweep up to determine number of proxy points and shells - do this with source points in a box half the size of the spreading box
* Sweep down to determine nspread
* Determine bin side that is an integer multiple of dx closest to half the box size
* Sweep up proxy points to bring the error back down to tolerance


### Desired demos
* Chunkie demo
* FMM3DBIE demo
* Vector valued demo
* Gradient demo
* Flexural/4th order










