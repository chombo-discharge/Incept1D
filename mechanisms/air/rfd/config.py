# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Configuration class for the relativistic feedback discharge model.

Implements the protocol described in :ref:`Chap:ConfigPy`.  Three parameters
can be set, all needed while the module body runs: the photon attenuation
lengths, which depend on the horizontal radius of the high-field region, and
the sea-level density the coefficients are relative to.  Unknown keys are
refused, since a configuration written for another mechanism would otherwise
be solved as this one; keys starting with ``_`` are comments
(``_comment``) or added by the loader (``_cfg_dir``) and are always allowed.
"""

from dataclasses import dataclass
from typing import Optional

from incept1d.constants import kB


@dataclass
class Config:
    """One case of the relativistic feedback discharge model."""

    label: str = ""
    #: Forward and backward photon attenuation lengths at sea level [m].
    photon_attenuation_fwd: Optional[float] = None
    photon_attenuation_bck: Optional[float] = None
    #: Pressure [bar] and temperature [K] of the sea-level reference density.
    sea_level_pressure: float = 1.01325
    sea_level_temperature: float = 293.15

    def pre_exec_vars(self):
        """The attenuation lengths and N0 are read at module level."""
        out = {"_N0": self.sea_level_pressure * 1e5 / (kB * self.sea_level_temperature)}
        if self.photon_attenuation_fwd is not None:
            out["_LAMBDA_A_FWD"] = float(self.photon_attenuation_fwd)
        if self.photon_attenuation_bck is not None:
            out["_LAMBDA_A_BCK"] = float(self.photon_attenuation_bck)
        return out

    def post_exec_init(self, mod):
        """Nothing to initialise."""

    def mechanism_params(self):
        """No run-time scalings."""
        return {}

    @classmethod
    def from_dict(cls, d, cfg_dir=None):
        """Build from a raw JSON dict; keys this mechanism does not know are errors."""
        known = {f.name for f in cls.__dataclass_fields__.values()}
        d = {k: v for k, v in d.items() if not k.startswith("_")}
        unknown = sorted(set(d) - known)
        if unknown:
            raise ValueError(
                f"Unknown key(s) {unknown} for the RFD mechanism; "
                f"expected some of {sorted(known)}"
            )
        return cls(**d)
