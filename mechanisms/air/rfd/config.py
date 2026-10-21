# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Configuration class for the relativistic feedback discharge model.

Implements the protocol described in :ref:`Chap:ConfigPy`.  Three parameters
can be set, all needed while the module body runs: the photon attenuation
lengths, which depend on the horizontal radius of the high-field region, and
the sea-level number density the relative density is taken against.  Unknown keys are
refused, since a configuration written for another mechanism would otherwise
be solved as this one; keys starting with ``_`` are comments
(``_comment``) or added by the loader (``_cfg_dir``) and are always allowed.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class Config:
    """One case of the relativistic feedback discharge model."""

    label: str = ""
    #: Forward and backward photon attenuation lengths at sea level [m].
    photon_attenuation_fwd: Optional[float] = None
    photon_attenuation_bck: Optional[float] = None
    #: Sea-level number density [m^-3], N0 of n~ = N / N0: the authors' value.
    sea_level_density: float = 2.688e25

    def pre_exec_vars(self):
        """The attenuation lengths and N0 are read at module level."""
        out = {"_N0": float(self.sea_level_density)}
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
