# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath("../../nsforest"))
sys.path.insert(0, os.path.abspath("../.."))

# -- Project information

project = 'NS-Forest'
copyright = '202X, JCVI, NIH'
author = 'Beverly Peng'

release = '4.1'
version = '4.1'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'nbsphinx', 
]

exclude_patterns = ['_build', '**.ipynb_checkpoints']

intersphinx_mapping = {
    'IPython': ('https://ipython.readthedocs.io/en/stable/', None),
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
    'anndata': ("https://anndata.readthedocs.io/en/stable/", None),
    'pandas': ("https://pandas.pydata.org/pandas-docs/stable/", None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
