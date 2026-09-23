# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Constants.py — physical constants used throughout the project.

All values are sourced from scipy.constants (CODATA 2018 recommended values).
"""

from scipy.constants import k as kB, e as Q, c as c_light

__all__ = ["kB", "Q", "c_light"]
