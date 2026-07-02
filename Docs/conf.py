"""Sphinx configuration for the Incept1D documentation."""

import os
import sys

# Make the repository root importable so autodoc can `import Inception`,
# `import Eigenvalues`, etc. exactly as they are invoked from the command
# line (no package/src layout is used in this project).
sys.path.insert(0, os.path.abspath(".."))

# Force a non-interactive matplotlib backend before autodoc imports any of
# the project modules (several of them do `import matplotlib.pyplot as plt`
# and call plotting functions at call-time; Agg avoids requiring a display
# in CI / pre-commit environments).
import matplotlib  # noqa: E402

matplotlib.use("Agg")

# -- Project information -----------------------------------------------------

project = "Incept1D"
copyright = "2026, Robert Marskar"
author = "Robert Marskar"
release = "0.1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Napoleon: the project uses NumPy-style docstrings throughout.
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
# Render numpydoc "Attributes" sections as a field list instead of separate
# `.. attribute::` directives, which would otherwise collide with the
# attribute directives autodoc already generates for annotated (e.g.
# dataclass) class attributes.
napoleon_use_ivar = True

# Autodoc defaults.
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

# Mock heavy/optional imports that are not needed to *read* the API but that
# may be unavailable in a minimal docs-build environment.
autodoc_mock_imports = []

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
