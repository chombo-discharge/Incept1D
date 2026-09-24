# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Configuration class for the toy mechanism.

Implements the same three-step protocol as ``mechanisms/air/air_config.py``, so the
test suite exercises the real ``incept1d.mechanism.load_mechanism`` path:

    config.pre_exec_vars()      -> dict injected before the module body runs
    config.post_exec_init(mod)  -> called after exec
    config.mechanism_params()   -> kwargs forwarded to Mechanism()
    config.label                -> display label
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Config:
    """Parameters defining one toy-mechanism configuration."""

    label: str = ""
    alpha_A: Optional[float] = None
    alpha_B: Optional[float] = None
    eta_0: Optional[float] = None
    delta_0: Optional[float] = None
    gamma0: Optional[float] = None
    xi_photo: float = 1.0
    xi_emit: float = 1.0
    reaction_multipliers: dict = field(default_factory=dict)
    pos_override: Optional[dict] = None
    neg_override: Optional[dict] = None
    # Set by post_exec_init so tests can assert that it ran.
    initialised: bool = False

    def pre_exec_vars(self):
        """Module-level names injected into the mechanism before it executes."""
        out = {}
        for attr, name in (
            ("alpha_A", "_ALPHA_A"),
            ("alpha_B", "_ALPHA_B"),
            ("eta_0", "_ETA_0"),
            ("delta_0", "_DELTA_0"),
            ("gamma0", "_GAMMA0"),
        ):
            value = getattr(self, attr)
            if value is not None:
                out[name] = value
        return out

    def post_exec_init(self, mod):
        """Hook run after the mechanism module body has executed."""
        self.initialised = True
        mod._POST_EXEC_RAN = True

    def mechanism_params(self):
        """Keyword arguments forwarded to ``Mechanism``."""
        return {
            "reaction_multipliers": dict(self.reaction_multipliers),
            "xi_photo": self.xi_photo,
            "xi_emit": self.xi_emit,
            "gamma0": self.gamma0,
            "pos_override": self.pos_override,
            "neg_override": self.neg_override,
        }

    @classmethod
    def from_dict(cls, d, cfg_dir=None):
        """Build a Config from a raw JSON-style dict, ignoring unknown keys."""
        known = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in known})
