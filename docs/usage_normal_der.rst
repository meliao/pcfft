Example usage 2
================

Suppose we want to evaluate the following sum:

.. math::    
    
   g(y_i) = \sum_{j=1}^N \partial_{\boldsymbol{n}_j} k(y_i, x_j) \mu_j

where :math:`k` is a 2D logarithmic kernel, which arises as a multiple of the Green's function of the Laplace equation in 2D:

.. math:: k(y_i-x_j) = \log ( \|y_i-x_j\| )

This type of sum arises in when solving a boundary integral equation formulation of an elliptic BVP with Dirichlet boundary conditions.
:doc:`usage` shows how to evaluate a similar sum without the normal derivatives. First, we have to define the kernel and its gradient. See the note in :doc:`api` for more details on how to specify kernels and points.

.. code:: matlab

   function k_evals = kern(src_pts, target_pts)
      % src_pts is a struct where src_pts.r with shape (2, M)
      % target_pts is a struct where target_pts.r has shape (2, N)
      % Computes log{|| src - target||}
      % Output shape is (N, M)

      % Shape (N, M)
      rx = src_pts.r(1, :) - target_pts.r(1, :).';
      ry = src_pts.r(2, :) - target_pts.r(2, :).';

      dist = sqrt(rx.^2 + ry.^2);

      k_evals = log(dist);
   end

   function k_evals = kern_s(src_pts, target_pts)
      % src_pts is a struct where 
      %         - src_pts.r with shape (2, M)
      %         - src_pts.n with shape (2, M)
      % target_pts is a struct where  target_pts.r has shape (2, N)
      % Computes \partial_{n(src)}log{|| src - target||}
      % Output shape is (N, M)

      % Shape (N, M)
      rx = src_pts.r(1, :) - target_pts.r(1, :).';
      ry = src_pts.r(2, :) - target_pts.r(2, :).';

      dist = sqrt(rx.^2 + ry.^2);

      k_evals = src_pts.n(1,:).*rx + src_pts.n(2,:).*ry;
      k_evals = k_evals ./ dist.^2;
   end

We also need to define the source and target points, and specify the normal vector at each source point:

.. code:: matlab

   N = 1000; % number of sources
   M = 3000; % number of targets
   src_info = struct();
   src_info.r = rand(2, N);
   src_info.n = randn(2, N); % normal vector at each source point
   targ_info = struct();
   targ_info.r = rand(2, M);

We construct the regular grid just as we did before -- we only need to specify the free-space kernel to get the grid right:

.. code:: matlab

   tol = 1e-6;
   [grid_info, proxy_info] = get_grid(@kern, src_info, targ_info, tol);

Because we are taking derivatives with respect to the source points, we need to specify the kernel and its gradient when we call :func:`get_spread` for the sources. The targets are not differentiated, so we can just specify the free-space kernel and don't need to pass a list of required point_info fields.

.. code:: matlab

   [A_spread_src, srt_info_src] = get_spread(@kern, @kern_s, src_info, ...
                                             grid_info, proxy_info, {'r','n'});
   [A_spread_targ, srt_info_targ] = get_spread(@kern, [], targ_info, ...
                                                grid_info, proxy_info);

The precorrected FFT algorithm approximates the field due to the dipoles in ``src_info`` by a collection of point charges on the regular grid. The field due to these point charges can be quickly computed using an FFT. We now precompute the Fourier transform of the free-space kernel:

.. code:: matlab

   kern_hat = get_kernhat(@kern, grid_info);

The approximation of each dipole in ``src_info`` by grid charges is only valid away from that dipole. 
However, the spreading matrices spread all dipoles to the grid, so immediately using the spreading matrices will approximate nearby interactions with large errors.
As the final precomputation step, we use :func:`get_addsub` to compute the corrections to fix the errors caused by ignoring this. 
When we call this function, we need to provide the free-space kernel and the kernel containing the true interaction between the all source and target (including any derivatives):

.. code:: matlab

   A_addsub = get_addsub(@kern, @kern_s, src_info, targ_info, grid_info, ...
                         proxy_info, srt_info_src, srt_info_targ, ...
                         A_spread_src, A_spread_targ);

Now that we have finished our precomputations, we can evaluate the sum by calling :func:`pcfft_apply`:

.. code:: matlab

   sigma = rand(N, 1); % source strengths
   u = pcfft_apply(sigma, A_spread_src, A_spread_targ, A_addsub, kern_hat);

And that's it! The `source repository <https://github.com/meliao/pcfft/tree/main/demos>`_ contains other examples, including differentiation of source and target points, and integration with `chunkIE <https://github.com/fastalgorithms/chunkie>`_.
