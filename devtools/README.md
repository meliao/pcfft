# pcfft/devtools

* `accuracy_tests/` contains tests for the accuracy of the entire PCFFT algorithm applied to 2D and 3D problems for a range of kernels. They can be thought of as integration tests.
* `convergence_plots/` contains scripts which run various parts of the PCFFT algorithm, and record + plot the convergence of various quantities like error, nshell, nproxy, and dx with the requested tolerance.
* `test/` are unit tests for the various helper functions in the PCFFT package.
* `timing_tests/` contain scripts that loop through problem sizes and compare the time required to solve Helmholtz BVPs using PCFFT, FLAM, and chunkIE+FMM. 
* `visual_checks/` are unit tests which don't assert anything, but are helpful for visually checking outputs of various parts of the code, i.e. grid point placement.