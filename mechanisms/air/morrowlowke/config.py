# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
config.py — configuration class for ``air_morrowlowke.py``.

The shared air configuration (``../air_config.py``), restricted to the keys
this mechanism understands.  Its transport is analytic, so there is no
``cross_sections`` key.
"""

import os

from incept1d.mechanism import load_helper

_air_config = load_helper(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "air_config.py")
)


class Config(_air_config.AirConfig):
    """One configuration of ``air_morrowlowke.py``."""

    KEYS = _air_config.COMMON_KEYS
