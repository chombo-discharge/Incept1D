"""Sphinx configuration for the Incept1D documentation.

Build from the repository root with

    python3 -m sphinx -b html Docs/source Docs/build/html

or ``make html`` from inside ``Docs/``.
"""

import os
import subprocess
import sys
import warnings

# autodoc imports the `incept1d` package.  Normally it is installed in the
# build environment (`pip install -e .`); putting `src/` first on sys.path
# lets the docs build from a bare checkout as well and guarantees that the
# documented code is the one in this working tree.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src"))

# Force a non-interactive matplotlib backend before autodoc imports any of
# the project modules (several of them do `import matplotlib.pyplot as plt`
# at import time; Agg avoids requiring a display in CI / pre-commit
# environments).
import matplotlib  # noqa: E402

matplotlib.use("Agg")

# A stale apt-installed ``mpl_toolkits`` next to a pip-installed matplotlib
# emits a harmless "Unable to import Axes3D" warning on every import; keep it
# out of the build log so real Sphinx warnings stay visible.
warnings.filterwarnings("ignore", message="Unable to import Axes3D")

# -- Project information -----------------------------------------------------

project = "Incept1D"
copyright = "2026, SINTEF Energy Research"
author = "Robert Marskar"
release = "0.1.0"

try:
    commit_id = (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=_REPO_ROOT)
        .strip()
        .decode("ascii")
    )
except Exception:  # not a git checkout (e.g. a source tarball)
    commit_id = "unknown"

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

master_doc = "index"

# Number figures/tables/equations so they can be referenced with :numref: /
# :eq: from any page.
numfig = True
math_numfig = True
numfig_secnum_depth = 2
math_eqref_format = "Eq. {number}"

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

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# MathJax macros for the notation used on the theory pages
# (\bm for bold matrices/vectors, \diag).
mathjax3_config = {
    "tex": {
        "macros": {
            "bm": ["\\boldsymbol{#1}", 1],
            "diag": "\\operatorname{diag}",
        }
    }
}

# Substitutions available on every page.
rst_epilog = f"""
.. |commit| replace:: ``{commit_id}``
.. |repo| replace:: https://github.com/chombo-discharge/Incept1D
"""

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["my_theme.css"]
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
}

# -- Options for LaTeX output ------------------------------------------------

latex_engine = "xelatex"
latex_elements = {
    "papersize": "a4paper",
    "pointsize": "10pt",
    "figure_align": "htb",
    # bm is used for bold math (\bm{R}, \bm{V}, ...) throughout the theory
    # pages.  xelatex + DejaVu fonts handle the Unicode (Greek letters,
    # arrows, superscripts) that the autodoc'd docstrings contain.
    "preamble": r"\usepackage{bm}\newcommand{\diag}{\operatorname{diag}}",
    "fontpkg": r"\setmainfont{DejaVu Serif}\setsansfont{DejaVu Sans}\setmonofont{DejaVu Sans Mono}",
}
latex_documents = [
    ("index", "Incept1D.tex", "Incept1D Documentation", author, "manual"),
]
