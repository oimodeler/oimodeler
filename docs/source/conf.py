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
import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
import pathlib
import sys
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.join(pathlib.Path(__file__).parents[2].resolve().as_posix(),"oimodeler"))
print("**************************************")
print(sys.path)
print("**************************************")
# -- Project information -----------------------------------------------------

project = 'oimodeler'
copyright = '2021, A. Meilland'
author = 'A. Meilland'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
 'sphinx.ext.duration',
 'sphinx.ext.doctest',
 'sphinx.ext.autodoc',
 'sphinx.ext.autosummary',
 'sphinx_rtd_theme',
]

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
html_logo = "../../images/logo.png"
html_favicon = "../../images/favicon.ico"
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'style_nav_header_background':'#eeeeee',
    'collapse_navigation': False,
}
html_css_files = ["custom.css"]
#numpydoc_class_members_toctree = False