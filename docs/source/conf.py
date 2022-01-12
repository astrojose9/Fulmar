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
import sys
sys.path.insert(0, os.path.abspath('/home/jrodrigues/Documents/PhD/fulmar'))

# -- Project information -----------------------------------------------------

project = 'FULMAR'
copyright = '2021, José Rodrigues'
author = 'José Rodrigues'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# Use Numpy-style docstrings and not Google-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'

# html_logo = "FULMAR_logo_name.png"
html_theme_options = {

    "light_logo": "FULMAR_logo_light_mode.png",
    "dark_logo": "FULMAR_logo_name.png",

    "sidebar_hide_name": True,
    "light_css_variables": {
        "color-brand-primary": "#02B01A",
        "color-brand-content": "#02B01A",
        "color-problematic": "#00b300",
    },
    "dark_css_variables": {

        "color-foreground-primary": "#ffffffcc",  # for main text and headings
        "color-foreground-secondary": "#9ca0a2",  # for secondary text
        "color-foreground-muted": "#818d86",  # for muted text
        "color-foreground-border": "#666666",  # for content borders

        #        "color-background-primary": "#2a2c2a",  # for content
        #        "color-background-secondary": "#222522",  # for navigation + ToC
        "color-background-primary": "#222522",  # for content
        "color-background-secondary": "#2a2c2a",  # for navigation + ToC
        "color-background-hover": "#1e2421ff",  # for navigation-item hover
        "color-background-hover--transparent": "#1e242100",
        "color-background-border": "#303335",  # for UI borders
        "color-brand-primary": "#02B01A",
        "color-brand-content": "#02B01A",
        "color-highlighted-background": "#086533",

        "color-guilabel-background": "#0863358",
        "color-guilabel-border": "#135f3980",

        "color-api-highlight-on-target": "#330033",
        "color-card-background": "#181a18",

        "color-problematic": "#51ee51",
    },
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['.._build/html/_static']


# intersphinx enables links to classes/functions in the packages defined here:
intersphinx_mapping = {'python': ('https://docs.python.org/3/', None),
                       'numpy': ('https://docs.scipy.org/doc/numpy/', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
                       'matplotlib': ('https://matplotlib.org', None),
                       'astropy': ('https://docs.astropy.org/en/latest/', None),
                       'lightkurve': ('https://docs.lightkurve.org/', None)}

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#
add_module_names = False
