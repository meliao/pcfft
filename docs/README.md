## Documentation

The documentation is available at [https://pcfft.readthedocs.io/en/latest/](https://pcfft.readthedocs.io/en/latest/).

## How to compile a local copy

First set up a python virtual environment and install the following packages:
```
pip install sphinx sphinxcontrib-matlabdomain sphinx_rtd_thetme sphinxcontrib-bibtex
```

Then to compile, run
```
cd docs
make html
open _build/html/index.html
```

