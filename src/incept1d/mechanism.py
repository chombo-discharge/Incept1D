# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Mechanism loading: the ``Mechanism`` wrapper, JSON configuration files and
``load_mechanism``.

A mechanism is a plain Python file (see ``mechanisms/air/air_pancheshnyi.py``)
that is executed by :func:`load_mechanism` and must expose the interface
listed in :data:`REQUIRED_ATTRS`.  A companion ``config.py`` in the same
directory (optional) implements the ``pre_exec_vars()`` /
``post_exec_init()`` / ``mechanism_params()`` protocol used to inject
configuration parameters from JSON files.
"""

import importlib.util
import json
import os


class Mechanism:
    """Loaded mechanism module with configuration baked in.

    Wraps a bare mechanism module and applies all configuration parameters
    (xi_photo, xi_emit, reaction_multipliers, gamma overrides) directly inside
    the interface methods.  Callers never need a separate modifier object.

    Call resolve(positive) to obtain a polarity-resolved copy with the
    appropriate cathode SEE overrides applied for that polarity.
    """

    def __init__(
        self,
        _mod,
        label="",
        reaction_multipliers=None,
        xi_photo=1.0,
        xi_emit=1.0,
        gamma0=None,
        gamma1=None,
        eref=None,
        beta=None,
        pos_override=None,
        neg_override=None,
    ):
        self._mod = _mod
        self.label = label
        self._reaction_multipliers = reaction_multipliers or {}
        self._xi_photo = float(xi_photo)
        self._xi_emit = float(xi_emit)
        self._gamma0 = gamma0
        self._gamma1 = gamma1
        self._eref = eref
        self._beta = beta
        self._pos_override = pos_override  # raw dict or None
        self._neg_override = neg_override  # raw dict or None
        self.SPECIES = _mod.SPECIES
        self.ELECTRON_INDEX = _mod.ELECTRON_INDEX
        if hasattr(_mod, "alpha"):
            self.alpha = lambda EN, p, T: _mod.alpha(EN, p, T)
        if hasattr(_mod, "eta"):
            self.eta = lambda EN, p, T: _mod.eta(EN, p, T)

    def resolve(self, positive: bool) -> "Mechanism":
        """Return a polarity-resolved copy with cathode SEE overrides applied."""
        ovr = self._pos_override if positive else self._neg_override
        if ovr is None:
            return self
        return Mechanism(
            self._mod,
            self.label,
            reaction_multipliers=self._reaction_multipliers,
            xi_photo=ovr.get("xi_photo", self._xi_photo),
            xi_emit=ovr.get("xi_emit", self._xi_emit),
            gamma0=ovr.get("gamma0", self._gamma0),
            gamma1=ovr.get("gamma1", self._gamma1),
            eref=ovr.get("eref", self._eref),
            beta=ovr.get("beta", self._beta),
        )

    def get_R(self, EN, p, T):
        return self._mod.get_R(EN, p, T, multipliers=self._reaction_multipliers)

    def get_V(self, EN, p, T):
        return self._mod.get_V(EN, p, T)

    def get_B(self, EN, p, T):
        return self._mod.get_B(EN, p, T) * self._xi_photo

    def get_C(self, EN, p, T):
        return self._mod.get_C(EN, p, T)

    def get_kappa(self, p, T):
        return self._mod.get_kappa(p, T)

    def get_gamma_plus(self, EN, p, T):
        return self._mod.get_gamma_plus_with(
            EN,
            p,
            T,
            gamma0=self._gamma0,
            gamma1=self._gamma1,
            eref=self._eref,
            beta=self._beta,
        )

    def get_gamma_Psi(self, EN, p, T):
        return self._mod.get_gamma_Psi(EN, p, T) * self._xi_emit

    def get_Pi_e(self):
        return self._mod.get_Pi_e()

    def get_Pi_plus(self):
        return self._mod.get_Pi_plus()

    def get_Pi_minus(self):
        return self._mod.get_Pi_minus()

    def init_photoionization(self, ngroups=3, cone_angle_deg=2.0):
        if hasattr(self._mod, "init_photoionization"):
            self._mod.init_photoionization(ngroups, cone_angle_deg)


#: Attributes that every mechanism module must expose.
REQUIRED_ATTRS = [
    "SPECIES",
    "ELECTRON_INDEX",
    "get_R",
    "get_V",
    "get_Pi_e",
    "get_Pi_plus",
    "get_Pi_minus",
    "get_gamma_plus",
    "get_gamma_plus_with",
    "get_B",
    "get_C",
    "get_kappa",
    "get_gamma_Psi",
]


def read_json_configs(json_paths):
    """Read one or more JSON config files and return a flat list of raw dicts.

    Each file may contain: a dict with a "configurations" key (list), a bare
    list of dicts, or a single dict.  Each returned dict has '_cfg_dir' set to
    the directory of the JSON file so that Config.from_dict can resolve
    relative paths (e.g. cross_sections) correctly.
    """
    results = []
    for path in json_paths:
        abs_path = os.path.abspath(path)
        cfg_dir = os.path.dirname(abs_path)
        with open(abs_path) as fh:
            data = json.load(fh)
        if isinstance(data, dict) and "configurations" in data:
            items = data["configurations"]
        elif isinstance(data, list):
            items = data
        else:
            items = [data]
        for item in items:
            d = dict(item)
            d.setdefault("_cfg_dir", cfg_dir)
            results.append(d)
    return results


def load_mechanism(path, config_dict=None):
    """
    Load a reaction mechanism module from a Python source file.

    The module must expose the standard interface defined in REQUIRED_ATTRS.

    Parameters
    ----------
    path : str
        Absolute or relative path to the mechanism .py file.
    config_dict : dict or None
        Raw configuration dict (e.g. parsed from JSON).  If given, a Config
        object is constructed from it by the companion config.py that lives
        alongside the mechanism file.  Pass None for a baseline configuration.

    Returns
    -------
    Mechanism
    """
    abs_path = os.path.abspath(path)
    mech_dir = os.path.dirname(abs_path)
    if not os.path.isfile(abs_path):
        raise FileNotFoundError(f"Mechanism file not found: {abs_path}")

    # Discover config.py in the same directory as the mechanism.
    cfg_py = os.path.join(mech_dir, "config.py")
    if os.path.isfile(cfg_py):
        _cspec = importlib.util.spec_from_file_location("_mechconfig", cfg_py)
        _cm = importlib.util.module_from_spec(_cspec)
        _cspec.loader.exec_module(_cm)
        config = _cm.Config.from_dict(config_dict or {}, mech_dir)
    elif config_dict:
        raise ImportError(
            f"No config.py found in {mech_dir} but a config_dict was supplied"
        )
    else:
        config = None

    spec = importlib.util.spec_from_file_location("mechanism", abs_path)
    mod = importlib.util.module_from_spec(spec)
    if config is not None:
        for attr, val in config.pre_exec_vars().items():
            mod.__dict__[attr] = val
    spec.loader.exec_module(mod)
    missing = [a for a in REQUIRED_ATTRS if not hasattr(mod, a)]
    if missing:
        raise AttributeError(
            f"Mechanism '{path}' is missing required attributes: {missing}"
        )
    if config is not None:
        config.post_exec_init(mod)
    label = config.label if config is not None else ""
    params = config.mechanism_params() if config is not None else {}
    return Mechanism(mod, label, **params)
