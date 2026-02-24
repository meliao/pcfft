.. pcfft documentation master file, created by
   sphinx-quickstart on Mon Feb  9 15:35:42 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pcfft documentation
===================

PCFFT is a MATLAB library for fast N-body summation of translation-invariant kernels in :math:`\mathbb{R}^2` and :math:`\mathbb{R}^3`. It can approximate sums of the form:

.. math:: 
   
   f(y_i ) = \sum_{j=1}^N k(y_i, x_j) \mu_j, \tag{1}

.. math::

   g(y_i) = \sum_{j=1}^N \partial_{\boldsymbol{n}_i} \partial_{\boldsymbol{n}_j} k(y_i, x_j) \mu_j. \tag{2}

at a large number of source points :math:`x_j` and target points :math:`y_i` efficiently. The code is designed to be easy to use and adaptable to a wide range of kernels, allowing the user to rapidly prototype large-scale numerical computations. In particular:

 * The method applies to a broad class of smooth translation-invariant kernels :math:`k(y, x) = k(y - x)`, not just kernels arising from the Green's function of an elliptic PDE. 
 * Analytical knowledge of the Fourier transform of the kernel is not required. 
 * After a one-time precomputation step, the apply step uses Fast Fourier Transforms and sparse linear algebra, making it very fast for large problems.


:doc:`usage` provides a brief introduction to using the package to evaluate sums of the form :math:`(1)`, and :doc:`usage_normal_der` shows how to evaluate sums of the form :math:`(2)`. More examples are being built at `<https://github.com/meliao/pcfft/tree/main/demos>`_.

Source repository
------------------
Available on GitHub at `<https://github.com/meliao/pcfft>`_.

Installation and Dependencies
------------------------------
PCFFT can be installed from source:

.. code:: bash

   git clone https://github.com/meliao/pcfft.git --recurse-submodules


The package depends on Kenneth Ho's package FLAM, available at `<https://github.com/fastalgorithms/FLAM>`_. 

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   usage
   usage_normal_der
   api
