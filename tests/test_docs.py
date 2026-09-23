# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Invariants of the Sphinx sources that a ``dummy`` build does not check.

The pre-commit hook and the CI ``rst`` job build with the ``dummy`` builder,
which resolves references but not images.  A figure path that is wrong only
in its number of ``../`` therefore passes both and fails the full HTML
build, which runs much later.  These tests close that gap without needing
Sphinx or the built figures, which are git-ignored build products.
"""

import re
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parent.parent / "docs"
SOURCE = DOCS / "source"

#: Every figure lands here; see ``docs/figures/Makefile``.
FIGURE_DIR = SOURCE / "figures"

_IMAGE_RE = re.compile(r"^\s*\.\.\s+(?:figure|image)::\s*(\S+)\s*$", re.M)


def _image_directives():
    """Yield (rst_path, target) for every figure and image directive."""
    for rst in sorted(SOURCE.rglob("*.rst")):
        for target in _IMAGE_RE.findall(rst.read_text()):
            yield rst, target


def _makefile_figures():
    """The figure stems ``docs/figures/Makefile`` knows how to build."""
    text = (DOCS / "figures" / "Makefile").read_text()
    stems = set()
    for var in ("FIGURES", "OPTIONAL_FIGURES"):
        m = re.search(rf"^{var}\s*:=\s*(.*)$", text, re.M)
        assert m, f"{var} not found in docs/figures/Makefile"
        stems.update(m.group(1).split())
    return stems


def test_there_are_image_directives_to_check():
    """Guard against the collection silently matching nothing."""
    assert list(_image_directives())


@pytest.mark.parametrize("rst,target", list(_image_directives()), ids=str)
def test_image_resolves_into_the_figure_directory(rst, target):
    """D1: a figure path must land in docs/source/figures/.

    This is what catches a stale ``../`` after a page is moved between
    chapters: the path still parses, but points at a directory that holds
    no figures and never will.
    """
    resolved = (rst.parent / target).resolve()
    assert resolved.parent == FIGURE_DIR, (
        f"{rst.relative_to(SOURCE)} references {target}, which resolves to "
        f"{resolved.parent}, not {FIGURE_DIR}"
    )


@pytest.mark.parametrize("rst,target", list(_image_directives()), ids=str)
def test_image_is_a_figure_the_makefile_builds(rst, target):
    """D2: the figure must be one docs/figures/Makefile produces.

    The figures themselves are build products and are not in the
    repository, so the Makefile's own list is the reference.
    """
    stem = Path(target).name.rsplit(".", 1)[0]
    assert stem in _makefile_figures(), (
        f"{rst.relative_to(SOURCE)} references figure '{stem}', which "
        f"docs/figures/Makefile does not build"
    )
