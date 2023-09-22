# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import os
import sys
from pathlib import Path


print("Python version")
print (sys.version)
sys.path.insert(0,str(Path(__file__).parent.parent.parent ))
sys.path.insert(0,str(Path(__file__).parent.parent.parent / "oimodeler"))
sys.path.insert(0,str(Path("sphinxext").resolve()))
print(sys.path)

import oimodeler as oim
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.

# -- Project information -----------------------------------------------------

project = 'oimodeler'
release = oim.__version__
copyright = '2023, A. Meilland'
author = 'A. Meilland'

print(f"OIMODELER v{release}")
# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'autodocsumm',
    'sphinx_rtd_theme',
]

napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '_templates']

autosummary_generate = True
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_logo = "../../images/oimodelerlogo_small.png"
html_favicon = "../../images/favicon.ico"
html_theme_options = {
    'logo_only': True,
    'display_version': True,
    'style_nav_header_background': '#eeeeee',
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
}

html_context = {
    'css_files': ['_static/custom.css'],
}
autodoc_member_order = 'bysource'
# numpydoc_class_members_toctree = False
"""
autodoc_default_options = {
    'members': 'var1, var2',
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__',
    'show-inheritance': False
}
"""

autodoc_default_options = {'autosummary': True,
                           'autosummary-no-nesting':False}


