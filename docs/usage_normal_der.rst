Example usage 2
================

Suppose we want to evaluate the following sum:

.. math::    
    
   g(y_i) = \sum_{j=1}^N \partial_{\boldsymbol{n}_i} \partial_{\boldsymbol{n}_j} k(y_i, x_j) \mu_j

where :math:`k` evaluates the kernel for 3D electrostatic interactions:

.. math:: k(y-x_j) = \frac{1}{4\pi \|y-x_j\|}

:doc:`usage` shows how to evaluate this sum without the normal derivatives, and shows how to define the kernel.

