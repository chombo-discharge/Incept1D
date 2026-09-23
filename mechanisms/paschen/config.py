# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Configuration class for the classical Townsend mechanism.

Implements the protocol described in :ref:`Chap:ConfigPy`.  Two parameters are
worth changing: which gas to use, and the cathode secondary-emission yield.
Both must be known before the module body runs -- the gas selects the Townsend
coefficients at module level -- so both go through ``pre_exec_vars``.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Config:
    """One gas / cathode combination of the classical Townsend model."""

    label: str = ""
    gas: str = "air"
    gamma0: Optional[float] = None
    mu_e_N: Optional[float] = None
    mu_ion_N: Optional[float] = None
    reaction_multipliers: dict = field(default_factory=dict)
    pos_override: Optional[dict] = None
    neg_override: Optional[dict] = None

    def pre_exec_vars(self):
        """The gas and the yield are needed while the module executes."""
        out = {"GAS": self.gas}
        for attr, name in (
            ("gamma0", "_GAMMA0"),
            ("mu_e_N", "_MU_E_N"),
            ("mu_ion_N", "_MU_ION_N"),
        ):
            value = getattr(self, attr)
            if value is not None:
                out[name] = value
        return out

    def post_exec_init(self, mod):
        """Nothing to initialise: there is no photoionization model here."""

    def mechanism_params(self):
        """Run-time scalings applied inside the Mechanism accessors."""
        return {
            "reaction_multipliers": dict(self.reaction_multipliers),
            "pos_override": self.pos_override,
            "neg_override": self.neg_override,
        }

    @classmethod
    def from_dict(cls, d, cfg_dir=None):
        """Build from a raw JSON dict, ignoring keys this mechanism does not use."""
        known = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in known})
