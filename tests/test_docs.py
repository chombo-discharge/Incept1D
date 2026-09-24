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
import shlex
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
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


EXAMPLES = ROOT / "examples"

#: One directory per worked example, each with a README and a docs page.
EXAMPLE_DIRS = sorted(
    d for d in EXAMPLES.iterdir() if d.is_dir() and not d.name.startswith((".", "_"))
)

_BASH_BLOCK_RE = re.compile(r"^```bash\n(.*?)^```", re.M | re.S)


def _normalise(text):
    """Join backslash continuations and collapse whitespace."""
    return " ".join(re.sub(r"\\\s*\n", " ", text).split())


def _main_command(example):
    """The README's main command: its first ``bash`` block running incept1d."""
    readme = (example / "README.md").read_text()
    for block in _BASH_BLOCK_RE.findall(readme):
        if block.lstrip().startswith("incept1d "):
            return _normalise(block)
    return None


def test_there_are_examples_to_check():
    """Guard against the collection silently matching nothing."""
    assert len(EXAMPLE_DIRS) >= 5


@pytest.mark.parametrize("example", EXAMPLE_DIRS, ids=lambda p: p.name)
def test_example_is_a_readme_not_a_script(example):
    """D3: an example is documented commands, not a script to run blindly.

    The reader should see the command and change it, so the README carries
    it; the figure Makefile carries its own copy for the build.
    """
    assert (example / "README.md").is_file(), f"{example.name} has no README.md"
    assert not list(example.glob("*.sh")), f"{example.name} still has a script"
    assert _main_command(example), (
        f"examples/{example.name}/README.md has no ```bash block starting "
        f"with 'incept1d'"
    )


@pytest.mark.parametrize("example", EXAMPLE_DIRS, ids=lambda p: p.name)
def test_main_command_is_valid(example):
    """D4: the README command parses, and names files that exist.

    Parsing uses the real ``incept1d`` parser without running anything, so
    a renamed or removed option breaks this test rather than the reader.
    """
    from incept1d.cli import build_parser

    argv = shlex.split(_main_command(example))
    assert argv[0] == "incept1d"
    args = build_parser().parse_args(argv[1:])
    for path in [args.mechanism, *getattr(args, "configs", [])]:
        assert (ROOT / path).is_file(), f"{example.name}: {path} does not exist"
    out = getattr(args, "write_to_file", None)
    assert out is None or Path(out).parent == Path("examples") / example.name, (
        f"{example.name}: the main command writes to {out}, not into "
        f"examples/{example.name}/"
    )


@pytest.mark.parametrize("example", EXAMPLE_DIRS, ids=lambda p: p.name)
def test_docs_page_shows_the_main_command(example):
    """D5: the documentation page runs exactly what the README runs."""
    page = SOURCE / "examples" / f"{example.name}.rst"
    assert page.is_file(), f"no docs page {page.relative_to(ROOT)}"
    assert _main_command(example) in _normalise(page.read_text()), (
        f"{page.relative_to(ROOT)} does not show the main command of "
        f"examples/{example.name}/README.md"
    )
