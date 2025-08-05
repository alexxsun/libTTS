# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'libtts'
copyright = '2025, Xin Xu'
author = 'Xin Xu'
release = '0.01'

# -- Configuration ---------------------------------------------------
import os
import sys
#sys.path.insert(0, os.path.abspath('../../src'))
# needed for editable installment
sys.path.insert(0, os.path.abspath('../../../build/cp312-cp312-linux_x86_64'))
import importlib.util

# Try to locate the installed libtts module
spec = importlib.util.find_spec("libtts")
if spec and spec.origin:
    libtts_dir = os.path.dirname(spec.origin)
    print("libtts_dir:", libtts_dir)
    sys.path.insert(0, libtts_dir)


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',  # Core extension to pull docs from docstrings
    'sphinx.ext.napoleon', # Teaches Sphinx to understand Google-style docstrings
    'sphinx.ext.viewcode', # Adds a '[source]' link to your function signatures
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
