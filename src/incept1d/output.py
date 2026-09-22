# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Shared helpers for the tab-separated result files written by the
``--write-to-file`` options of the command-line tools.

Every result file starts with a metadata block recording when, from which
git revision and with which command line it was produced, so that a data
file can always be traced back to the code that made it.
"""

import datetime
import os
import subprocess
import sys

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))


def git_revision(extra_paths=()):
    """
    Return a short description of the git revision of the installed package.

    Parameters
    ----------
    extra_paths : iterable of str
        Additional files (e.g. the mechanism file) whose uncommitted
        modifications should also mark the revision as dirty.

    Returns
    -------
    str
        ``"<short hash>"``, ``"<short hash> (dirty)"`` if the package sources
        or any of *extra_paths* have uncommitted changes, or
        ``"unavailable"`` when the package does not run from a git checkout.
    """
    try:
        git_hash = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=_PKG_DIR,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        paths = [_PKG_DIR] + [os.path.abspath(p) for p in extra_paths]
        dirty = (
            subprocess.check_output(
                ["git", "status", "--porcelain", "--"] + paths,
                cwd=_PKG_DIR,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        return git_hash + (" (dirty)" if dirty else "")
    except Exception:
        return "unavailable"


def write_metadata_header(fh, extra_paths=()):
    """
    Write the common ``# --- METADATA ---`` block to an open text file.

    Records the timestamp, the git revision (see :func:`git_revision`) and
    the full command line.  Callers append their own tool-specific header
    lines afterwards.
    """
    fh.write("# --- METADATA ---\n")
    fh.write(f"# Date:    {datetime.datetime.now().isoformat(timespec='seconds')}\n")
    fh.write(f"# Git:     {git_revision(extra_paths)}\n")
    fh.write(f"# Command: {' '.join(sys.argv)}\n")
    fh.write("# ---\n")
