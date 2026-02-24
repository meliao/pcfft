# pcfft
This repository provides routines for quickly computing N-Body calculations in 2D and 3D and is designed to work with generic translation invariant kernels

$$u_i = \sum_{j\neq i} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j.$$

The routines can also be used to compute sums involving source and target derivatives, such as

$$w_i = \sum_{j\neq i} \partial_{\mathbf{n}_i}\partial_{\mathbf{n}_j} K(\mathbf{x}_i - \mathbf{y}_j)\mu_j.$$

**NOTE:** The routines always skip terms with $x_i=y_j$.

## Installation instructions

PCFFT can be installed from source
```
git clone https://github.com/meliao/pcfft.git --recurse-submodules
```

The `--recurse-submodules` flag ensures the [FLAM](https://github.com/fastalgorithms/FLAM) package is available.

## User-callable routines

Detailed documentation is being built in `docs/`.

```
[grid_info, proxy_info] = get_grid(kernel, src_info, targ_info, ...
        tol, n_nbr)
```
This routine determines the size of the equispaced grid, the number of proxy points used for spreading, and the spreading parameters.

```
[A_spread, sort_info] = get_spread(kern_0, kern_der, ...
                                            src_info, grid_info, proxy_info, der_fields)
```
This routine returns the matrix that maps charge strengths at `src_info.r` to charge strengths on the equispaced grid, it also returns some point binning info used in `get_addsub()`

```
A_addsub = get_addsub(kern_0, kern_st, src_info, ...
    targ_info, grid_info, proxy_info, sort_info_s, sort_info_t, ...
    A_spread_s, A_spread_t)
```
This routine returns the matrix that fixes the near-field interactions that are done incorrectly by the spreading.

```
kern_hat = get_kernhat(kern_0, grid_info)
```
Evaluate FFT of kern_0 on the equispaced grid.

```
u = pcfft_apply(sigma, A_spread_s, A_spread_t, A_addsub, kern_0hat)
```
Compute the N-body sum using a precorrected FFT.

## Sample use

See the `demos/` directory for a variety of examples, including demonstrating pairing this package with some popular repositories.


## wishlist:
* vector valued kernels <- after profiling
* Profiling, particularly on an old computer


## proposed parameter plan:

* Determine spreading box half side length
* Sweep up to determine number of proxy points and shells - do this with source points in a box half the size of the spreading box
* Sweep down to determine nspread
* Determine bin side that is an integer multiple of dx closest to half the box size
* Sweep up proxy points to bring the error back down to tolerance


## Desired demos
* Chunkie demo -- separate demos for basic use+plotting and time comparison with FMM.
* FMM3DBIE demo
* Vector valued demo
* Gradient demo
* Flexural/4th order
* Matern kernel



## Compiling docs

First set up a python virtual environment and install the following packages:
```
pip install sphinx sphinxcontrib-matlabdomain sphinx_rtd_thetme
```

Then to compile, do
```
cd docs
make html
open _build/html/index.html
```





