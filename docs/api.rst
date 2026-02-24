PCFFT API
===================

Specifying kernels and points
------------------------------

.. note:: @Tristan can you go over this section and fill in what is missing? Things like `src_info.n`, kernel returning grad, ... ?

To specify points, the user should construct a structure with field ``r`` containing the coordinates of the points. If there are :math:`N` source points in :math:`\mathbb{R}^3`, the ``r`` field should be an array with shape ``[3, N]``. The same applies to target points.

.. code:: matlab

   N = 1000; % number of sources
   M = 3000; % number of targets
   src_info = struct();
   targ_info = struct();
   src_info.r = rand(3, N);
   targ_info.r = rand(3, M);


The PCFFT package requires the kernel to be specified as a MATLAB function handle. We assume a few things about the kernel function:

 * The first argument is the source point structure, and the second argument is the target point structure. The kernel function should not require other arguments.
 * If there are :math:`N` source points and :math:`M` target points, the kernel function should return an array of shape ``[M, N]`` containing the kernel evaluations between each target and source point.

Here is an example of a kernel function which matches these requirements:

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

In the following documentation, we will use the types ``point_info`` and  ``kernel`` to refer to structures containing point information and kernel function handles, respectively.

Function API
-------------
.. mat:currentmodule:: pcfft
.. mat:autofunction:: get_grid
.. mat:autofunction:: get_spread
.. mat:autofunction:: get_addsub
.. mat:autofunction:: get_kernhat
.. mat:autofunction:: pcfft_apply


Built-in classes
-------------------
.. autoclass:: pcfft.classes.GridInfo
   :members:

.. autoclass:: pcfft.classes.ProxyInfo
   :members:

.. autoclass:: pcfft.classes.SortInfo
   :members: