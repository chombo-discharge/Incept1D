# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
config.py — configuration class for Air mechanism files.

Defines Config for use with the air mechanism files in this directory.
Implements the config protocol expected by incept1d.mechanism.load_mechanism:

    config.pre_exec_vars()      → dict of {attr: value} to inject before exec
    config.post_exec_init(mod)  → called after exec (e.g. init_photoionization)
    config.mechanism_params()   → dict of kwargs forwarded to Mechanism()
    config.label                → str display label
"""

import os
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Config:
    """All parameters that define one Air-model configuration."""

    label: str = "Baseline"
    cross_sections: Optional[str] = None
    ngroups: int = 3
    cone_angle: float = 45.0
    xi_photo: float = 1.0
    xi_emit: float = 1.0
    reaction_multipliers: dict = field(default_factory=dict)
    gamma0: Optional[float] = None
    gamma1: Optional[float] = None
    eref: Optional[float] = None
    beta: Optional[float] = None
    positive_override: Optional[dict] = None
    negative_override: Optional[dict] = None

    # ── config protocol ──────────────────────────────────────────────────────

    def pre_exec_vars(self) -> dict:
        """Variables to inject into the module namespace before exec."""
        d = {}
        if self.cross_sections is not None:
            d["BOLSIG_FILE"] = self.cross_sections
        for attr, val in [
            ("_GAMMA0", self.gamma0),
            ("_GAMMA1", self.gamma1),
            ("_EREF", self.eref),
            ("_BETA", self.beta),
        ]:
            if val is not None:
                d[attr] = val
        return d

    def post_exec_init(self, mod) -> None:
        """Post-exec initialisation (photoionization group setup)."""
        if hasattr(mod, "init_photoionization"):
            mod.init_photoionization(self.ngroups, self.cone_angle)

    def mechanism_params(self) -> dict:
        """Kwargs forwarded to Mechanism() for runtime scaling."""
        return dict(
            reaction_multipliers=self.reaction_multipliers,
            xi_photo=self.xi_photo,
            xi_emit=self.xi_emit,
            pos_override=self.positive_override,
            neg_override=self.negative_override,
        )

    # ── construction from raw dict ────────────────────────────────────────────

    @classmethod
    def from_dict(cls, d: dict, mech_dir: str) -> "Config":
        """Construct Config from a raw JSON dict.

        Relative cross_sections paths are resolved first against
        d['_cfg_dir'] (injected by _read_json_configs), then mech_dir.
        """
        cs = d.get("cross_sections")
        if cs is not None and not os.path.isabs(cs):
            cfg_dir = d.get("_cfg_dir", mech_dir)
            candidate = os.path.join(cfg_dir, cs)
            if not os.path.isfile(candidate):
                candidate = os.path.join(mech_dir, cs)
            cs = candidate
        return cls(
            label=d.get("label", "Baseline"),
            cross_sections=cs,
            ngroups=d.get("ngroups", 3),
            cone_angle=d.get("cone_angle", 45.0),
            xi_photo=d.get("xi_photo", 1.0),
            xi_emit=d.get("xi_emit", 1.0),
            reaction_multipliers=d.get("reaction_multipliers", {}),
            gamma0=d.get("gamma0"),
            gamma1=d.get("gamma1"),
            eref=d.get("eref"),
            beta=d.get("beta"),
            positive_override=d.get("positive"),
            negative_override=d.get("negative"),
        )
