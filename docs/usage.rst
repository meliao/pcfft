Example usage 1
================

First, suppose we want to evaluate the following sum:

.. math:: f(y) = \sum_{j=1}^N k(y - x_j) \mu_j

where :math:`k` evaluates the kernel for 3D electrostatic interactions:

.. math:: k(y-x_j) = \frac{1}{4\pi \|y-x_j\|}

In this setting, we can think of the points :math:`x_j` as sources, each with strength :math:`\mu_j`. The points :math:`y_i` are targets. The goal is to evaluate :math:`f(y_i)` for many target points :math:`y_i` efficiently. We will call the number of target points :math:`M` and the number of source points :math:`N`.

First, we need to define the kernel as a function handle:

.. code:: matlab

   function evals = kern(src_info, targ_info)
       % evaluate electrostatic kernel between N source points src_info.r and 
       % M target points targ_info.r
       %
       % Output shape is (M, N)

       % Shape (M, N)
       rx = src_info.r(1, :) - targ_info.r(1, :).';
       ry = src_info.r(2, :) - targ_info.r(2, :).';
       rz = src_info.r(3, :) - targ_info.r(3, :).';

       dist = sqrt(rx.^2 + ry.^2 + rz.^2);
       k_evals = 1 ./ dist;
       evals = k_evals / (4*pi);
    end

Next, we need to define the source and target points. For simplicity, we will just use random points in a box:

.. code:: matlab

   N = 1000; % number of sources
   M = 3000; % number of targets
   src_info.r = rand(3, N);
   targ_info.r = rand(3, M);
   mu = rand(N, 1);

To evaluate the sum, we first need to construct a regular grid that will be used for the FFT. We need to choose the grid spacing carefully so that sources spread on the regular grid can accurately represent the far-field kernel evaluations. :func:`get_grid` handles this for us:

.. code:: matlab

   tol = 1e-6;
   [grid_info, proxy_info] = get_grid(kern, src_info, targ_info, tol);

Next, we call :func:`get_spread` to compute the matrices which spread the sources and targets to the regular grid.

.. code:: matlab

   [A_spread_src, ~, srt_info_src] = get_spread(kern, [], src_info, ...
                                                grid_info, proxy_info);
   [A_spread_targ, ~, srt_info_targ] = get_spread(kern, [], targ_info, ...
                                                grid_info, proxy_info);

The matrices ``A_spread_src`` and ``A_spread_targ`` take care of the far-field interactions, but we need to correct for near-field interactions which must be computed exactly. We do this by calling :func:`get_addsub`:

.. code:: matlab

   A_addsub = get_addsub(kern, [], src_info, targ_info, grid_info, ...
                         proxy_info, srt_info_src, srt_info_targ, ...
                         A_spread_src, A_spread_targ);

Finally, we can evaluate the sum by calling :func:`get_kernhat` and :func:`pcfft_apply`:

.. code:: matlab

   kern_hat = get_kernhat(kern, grid_info);
   f = pcfft_apply(mu, A_spread_src, A_spread_targ, ...
                   A_addsub, kern_hat);

Notice that we didn't use ``mu`` until the final line of this example. If the source points, target points, and kernel are fixed, we can think of all of the steps before the call to :func:`pcfft_apply` as a precomputation step. Once the precomputation is done, we can apply the same operator to different source strengths :math:`\mu` very quickly by calling :func:`pcfft_apply` again with the new source strengths and the same matrices.