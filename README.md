# pcfft
This repository provides routines for quickly computing N-Body calculations in 2D and 3D and is designed to work with generic translation invariant kernels

$$u_i = \sum_{j\neq i} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j$$

Also supports sums of the form

$$v_i = \sum_{j\neq i} \partial_{\mathbf{n}_j} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j$$

$$w_i = \sum_{j\neq i} \partial_{\mathbf{n}_i} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j$$

$$\tilde{w}_i = \sum_{j\neq i} \partial_{\mathbf{n}_i}\partial_{\mathbf{n}_j} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j$$

### User-callable precomputation routines
[grid_info, pxyinfo] = get_grid(kern, srcinfo, targinfo, eps, fill)
This routine determines the size of the equispaced grid, the number of proxy points used for spreading, and the spreading parameters

[A_spread, sortinfo] = get_spread(kern_0, kern, srcinfo, grid_info, pxyinfo)
This routine returns the matrix that maps charge strengths at srcinfo.r to charge strengths on the equispaced grid, it also returns some point binning info used in get_addsub

A_add_sub = get_addsub(kern_0, kern_s, kern_t, kern_st, srcinfo, targinfo, grid_info, pxyinfo, sortinfo_s, sortinfo_t)
This routine returns the matrix that fixes the interactions that are done incorrectly by the fast apple

kern_hat = get_kernhat(kern_0, grid_info)
evaluate FFT of kern_0 on grid


### Sample use



```
str = A_spread_s * mu; str = full(str);

str_hat = fft2(reshape(str,size(kern_hat));
u_hat = kern_hat .* str_hat;
ugrid = ifft2(u_hat);
    
u = A_spread_t.'*ugrid(:) + A_add_sub*mu;
```



### Depenencies

 * [FLAM](https://github.com/klho/FLAM)


### wishlist:
vector valued kernels
kernels with source and target information
concentric proxy surfaces







