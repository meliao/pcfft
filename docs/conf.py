# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pcfft'
copyright = '2026, Tristan Goodwill, Owen Melia'
author = 'Tristan Goodwill, Owen Melia'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

extensions = ['sphinx.ext.autodoc', 'sphinxcontrib.matlab']
primary_domain = "mat"


# Options for sphinxcontrib-matlab
matlab_src_dir = '..'
autoclass_content = 'both'
matlab_short_links = True
