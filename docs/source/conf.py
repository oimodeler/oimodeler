# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
from pathlib import Path

print("Python version")
print (sys.version)
sys.path.insert(0,str(Path(__file__).parent.parent.parent ))
sys.path.insert(0,str(Path(__file__).parent.parent.parent / "oimodeler"))
sys.path.insert(0,str(Path(__file__).parent.parent.parent / "examples" / "notebooks"))
sys.path.insert(0,str(Path("sphinxext").resolve()))
print(sys.path)


project = 'oimodeler'
copyright = '2024, Anthony Meilland'
author = 'Anthony Meilland'
release = 'v0.9x'


extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'autodocsumm',
    'sphinx_rtd_theme',
    'matplotlib.sphinxext.plot_directive',
    #'myst_parser',
    #'sphinx_gallery.gen_gallery',
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


autodoc_default_options = {'autosummary': True,
                           'autosummary-no-nesting':False}

#sphinx_gallery_conf = {
#    'examples_dirs': ['../../examples/tests'],
#    'gallery_dirs': ['auto_examples'],
#}