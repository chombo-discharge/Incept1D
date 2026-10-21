# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Declarative mechanism files: a JSON description compiled into the solver
interface.

A mechanism file names the species, their charges and mobilities, the
reactions with their rate coefficients, the Townsend coefficients and the
electrode yields, with every coefficient a formula
(:class:`incept1d.expressions.Expression`).  :func:`load_json_mechanism`
checks the file strictly — an unknown key, an unknown name in a formula or a
reaction that does not balance is an error, not a silent default — and
returns a :class:`DeclarativeMechanism`, which has the interface the solver
uses (see :data:`incept1d.mechanism.REQUIRED_ATTRS`).  What the Python
interface asks a mechanism to build by hand follows from the file: the row
selectors from the charges, the signs of the drift velocities from the
charges, the reaction matrix from the reaction strings.

A configuration (``read_json_configs``) sets the file's parameters by name,
scales reactions by their reaction string (``reaction_multipliers``) and
overrides parameters for one polarity (``positive``, ``negative``).

Formulas are evaluated with these variables:

    EN   reduced field [Td]          E  field [V/m]
    N    neutral number density [m^-3]
    p    pressure [bar]              T  gas temperature [K]
    kB   Boltzmann constant [J/K]    qe elementary charge [C]
    x_<X>  molar fraction of gas component X
    every parameter and every definition, in file order
    alpha, eta (when the file gives them), in rates and yields
"""

import json
import os
import re

import numpy as np

from incept1d.constants import Q, kB
from incept1d.expressions import FUNCTIONS, Expression
from incept1d.reactions import _normalize, parse_reaction

#: The ``format`` a mechanism file must declare.
FORMAT = "incept1d-mechanism/1"

_TOP_KEYS = {
    "format",
    "description",
    "parameters",
    "gas",
    "untracked",
    "definitions",
    "species",
    "alpha",
    "eta",
    "reactions",
    "emission",
}
_SPECIES_KEYS = {"name", "charge", "mobility"}
_REACTION_KEYS = {"reaction", "rate"}
_EMISSION_KEYS = {"ions"}
_CONFIG_KEYS = {"label", "reaction_multipliers", "positive", "negative"}
_BUILTINS = ("EN", "E", "N", "p", "T", "kB", "qe")
_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")

#: The tracked species that is the electron.
ELECTRON = "e"
#: The name that stands for the total neutral density in a reaction.
TOTAL_DENSITY = "M"


class MechanismFileError(ValueError):
    """A mechanism file or configuration that does not describe a mechanism."""


def is_mechanism_file(path):
    """True if *path* is a JSON file declaring the mechanism :data:`FORMAT`."""
    if not path.endswith(".json"):
        return False
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return False
    return isinstance(data, dict) and data.get("format") == FORMAT


def load_json_mechanism(path, config_dict=None):
    """
    Read, check and compile a mechanism file.

    Parameters
    ----------
    path : str
        The ``.json`` mechanism file.
    config_dict : dict or None
        One configuration (as returned by
        :func:`incept1d.mechanism.read_json_configs`), or None for the
        file's own parameters.

    Returns
    -------
    DeclarativeMechanism
    """
    with open(path) as fh:
        try:
            spec = json.load(fh)
        except ValueError as exc:
            raise MechanismFileError(f"{path}: not valid JSON: {exc}") from None
    return DeclarativeMechanism(spec, config_dict, source=path)


def _keys(obj, allowed, where):
    """Check that *obj* is a dict whose keys are *allowed* (or start with _)."""
    if not isinstance(obj, dict):
        raise MechanismFileError(f"{where}: expected an object, got {obj!r}")
    unknown = sorted(k for k in obj if k not in allowed and not k.startswith("_"))
    if unknown:
        raise MechanismFileError(
            f"{where}: unknown key(s) {unknown}; expected some of {sorted(allowed)}"
        )


def _number(value, where):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise MechanismFileError(f"{where}: expected a number, got {value!r}")
    return float(value)


class DeclarativeMechanism:
    """
    A mechanism compiled from a declarative description.

    Parameters
    ----------
    spec : dict
        The parsed mechanism file.
    config : dict or None
        One configuration: parameter values, ``label``,
        ``reaction_multipliers``, ``positive`` / ``negative``.
    source : str
        Where *spec* came from, for error messages.
    """

    def __init__(self, spec, config=None, source="<mechanism>"):
        self._spec = spec
        self._source = source
        _keys(spec, _TOP_KEYS, source)
        if spec.get("format") != FORMAT:
            raise MechanismFileError(
                f"{source}: 'format' must be {FORMAT!r}, got {spec.get('format')!r}"
            )
        self._parameters = self._read_parameters(spec.get("parameters", {}))
        config = dict(config or {})
        _keys(config, _CONFIG_KEYS | set(self._parameters), f"{source}: configuration")
        self.label = str(config.get("label", ""))
        self._values = dict(self._parameters)
        self._values.update(self._overrides(config, "configuration"))
        self._polarity = {}
        for side in ("positive", "negative"):
            if config.get(side):
                where = f"{source}: '{side}' override"
                _keys(config[side], set(self._parameters), where)
                self._polarity[side] = self._overrides(config[side], where)
        self._compile()
        self._multipliers = self._read_multipliers(config.get("reaction_multipliers"))
        self._last = None

    # --- reading -------------------------------------------------------------

    def _read_parameters(self, params):
        if not isinstance(params, dict):
            raise MechanismFileError(f"{self._source}: 'parameters' must be an object")
        out = {}
        for name, value in params.items():
            if name.startswith("_"):
                continue
            self._check_name(name, "parameter")
            out[name] = _number(value, f"{self._source}: parameter {name!r}")
        return out

    def _overrides(self, values, where):
        """The parameter values in *values* (whose keys are already checked)."""
        return {
            k: _number(v, f"{where} {k!r}")
            for k, v in values.items()
            if k in self._parameters
        }

    def _check_name(self, name, what):
        if not _NAME.match(name) or name in _BUILTINS or name in FUNCTIONS:
            raise MechanismFileError(
                f"{self._source}: {what} name {name!r} must start with a letter, "
                f"contain only letters, digits and _, and not be one of "
                f"{sorted(set(_BUILTINS) | set(FUNCTIONS))}"
            )

    def _read_multipliers(self, multipliers):
        if multipliers is None:
            return {}
        if not isinstance(multipliers, dict):
            raise MechanismFileError(
                f"{self._source}: 'reaction_multipliers' must be an object"
            )
        known = {_normalize(r["reaction"]): i for i, r in enumerate(self._reactions)}
        out, unknown = {}, []
        for key, value in multipliers.items():
            index = known.get(_normalize(key))
            if index is None:
                unknown.append(key)
            else:
                out[index] = _number(value, f"{self._source}: multiplier {key!r}")
        if unknown:
            raise MechanismFileError(
                f"Configuration '{self.label}' sets multipliers for reactions that "
                f"'{os.path.basename(self._source)}' does not define: {unknown}"
            )
        return out

    # --- compiling -----------------------------------------------------------

    def _compile(self):
        spec, src = self._spec, self._source
        names = list(_BUILTINS) + list(self._parameters)

        gas = spec.get("gas", {})
        if not isinstance(gas, dict):
            raise MechanismFileError(f"{src}: 'gas' must be an object")
        self._gas = {}
        for comp, frac in gas.items():
            if comp.startswith("_"):
                continue
            self._gas[comp] = Expression(
                frac, list(self._parameters), f"{src}: gas {comp!r}"
            )
        fractions = {c: e(self._values) for c, e in self._gas.items()}
        if fractions and abs(sum(fractions.values()) - 1.0) > 1e-6:
            raise MechanismFileError(
                f"{src}: gas fractions {fractions} do not sum to 1"
            )
        self._fraction_names = {c: f"x_{c}" for c in self._gas}
        for c, n in self._fraction_names.items():
            self._check_name(n, f"gas fraction of {c!r}")
        names += list(self._fraction_names.values())

        self._definitions = []
        for name, text in (spec.get("definitions") or {}).items():
            if name.startswith("_"):
                continue
            self._check_name(name, "definition")
            if name in names:
                raise MechanismFileError(
                    f"{src}: definition {name!r} is already defined"
                )
            self._definitions.append(
                (name, Expression(text, names, f"{src}: definition {name!r}"))
            )
            names.append(name)

        self._alpha = self._eta = None
        if "alpha" in spec:
            self._alpha = Expression(spec["alpha"], names, f"{src}: alpha")
        if "eta" in spec:
            self._eta = Expression(spec["eta"], names, f"{src}: eta")
        names += [
            n for n, e in (("alpha", self._alpha), ("eta", self._eta)) if e is not None
        ]
        if self._alpha is not None:
            self.alpha = lambda EN, p, T: self._env(EN, p, T)["alpha"]
        if self._eta is not None:
            self.eta = lambda EN, p, T: self._env(EN, p, T)["eta"]

        species = spec.get("species")
        if not isinstance(species, list) or not species:
            raise MechanismFileError(f"{src}: 'species' must be a non-empty list")
        self.SPECIES, self._charge, self._mobility = [], [], []
        for i, sp in enumerate(species):
            where = f"{src}: species[{i}]"
            _keys(sp, _SPECIES_KEYS, where)
            for k in _SPECIES_KEYS:
                if k not in sp:
                    raise MechanismFileError(f"{where}: missing {k!r}")
            name = sp["name"]
            if not isinstance(name, str) or not name or name in self.SPECIES:
                raise MechanismFileError(f"{where}: name {name!r} is empty or repeated")
            charge = sp["charge"]
            if isinstance(charge, bool) or not isinstance(charge, int) or charge == 0:
                raise MechanismFileError(
                    f"{where}: charge of {name!r} must be a non-zero integer: "
                    f"a neutral species does not drift"
                )
            self.SPECIES.append(name)
            self._charge.append(charge)
            self._mobility.append(
                Expression(sp["mobility"], names, f"{where}: mobility")
            )
        if (
            ELECTRON not in self.SPECIES
            or self._charge[self.SPECIES.index(ELECTRON)] != -1
        ):
            raise MechanismFileError(
                f"{src}: the electron {ELECTRON!r} with charge -1 is required"
            )
        self.ELECTRON_INDEX = self.SPECIES.index(ELECTRON)
        self._sign = -np.sign(np.array(self._charge, dtype=float))
        self._positive = [i for i, z in enumerate(self._charge) if z > 0]
        self._negative = [
            i for i, z in enumerate(self._charge) if z < 0 and i != self.ELECTRON_INDEX
        ]

        untracked = spec.get("untracked", [])
        if not isinstance(untracked, list) or not all(
            isinstance(u, str) for u in untracked
        ):
            raise MechanismFileError(f"{src}: 'untracked' must be a list of names")
        self._reactions = []
        for i, rx in enumerate(spec.get("reactions") or []):
            where = f"{src}: reactions[{i}]"
            _keys(rx, _REACTION_KEYS, where)
            for k in _REACTION_KEYS:
                if k not in rx:
                    raise MechanismFileError(f"{where}: missing {k!r}")
            self._reactions.append(
                self._compile_reaction(rx, set(untracked), names, where)
            )

        emission = spec.get("emission", {})
        _keys(emission, _EMISSION_KEYS, f"{src}: emission")
        ions = emission.get("ions", {})
        _keys(ions, {self.SPECIES[i] for i in self._positive}, f"{src}: emission.ions")
        missing = [
            self.SPECIES[i] for i in self._positive if self.SPECIES[i] not in ions
        ]
        if missing:
            raise MechanismFileError(
                f"{src}: emission.ions gives no yield for {missing}; write 0 for none"
            )
        self._gamma_plus = [
            Expression(
                ions[self.SPECIES[i]],
                names,
                f"{src}: emission.ions[{self.SPECIES[i]!r}]",
            )
            for i in self._positive
        ]

    def _compile_reaction(self, rx, untracked, names, where):
        text = rx["reaction"]
        if not isinstance(text, str) or "->" not in text:
            raise MechanismFileError(
                f"{where}: {text!r} is not a reaction 'A + B -> C'"
            )
        reactants, products = parse_reaction(text)
        tracked = [s for s in reactants if s in self.SPECIES]
        if len(tracked) != 1 or reactants[tracked[0]] != 1:
            raise MechanismFileError(
                f"{where}: {text!r} must have exactly one tracked reactant, once; "
                f"the others are gas components or {TOTAL_DENSITY!r}"
            )
        neutral = {}
        for name, count in reactants.items():
            if name in self.SPECIES:
                continue
            if name != TOTAL_DENSITY and name not in self._gas:
                raise MechanismFileError(
                    f"{where}: reactant {name!r} in {text!r} is neither a species, "
                    f"a gas component nor {TOTAL_DENSITY!r}"
                )
            neutral[name] = count
        for name in products:
            if (
                name not in self.SPECIES
                and name not in self._gas
                and name != TOTAL_DENSITY
                and name not in untracked
            ):
                raise MechanismFileError(
                    f"{where}: product {name!r} in {text!r} is not a species, a gas "
                    f"component or {TOTAL_DENSITY!r}; list it in 'untracked' if it "
                    f"is not followed"
                )
        charge = dict(zip(self.SPECIES, self._charge))
        q_in = sum(charge.get(s, 0) * c for s, c in reactants.items())
        q_out = sum(charge.get(s, 0) * c for s, c in products.items())
        if q_in != q_out:
            raise MechanismFileError(f"{where}: {text!r} does not conserve charge")
        delta = np.array(
            [products.get(s, 0) - reactants.get(s, 0) for s in self.SPECIES],
            dtype=float,
        )
        return {
            "reaction": text,
            "driver": self.SPECIES.index(tracked[0]),
            "delta": delta,
            "neutral": neutral,
            "rate": Expression(rx["rate"], names, f"{where}: rate"),
        }

    # --- evaluating ----------------------------------------------------------

    def _env(self, EN, p, T):
        """All variables at (EN, p, T); the last one is kept."""
        key = (EN, p, T)
        if self._last is not None and self._last[0] == key:
            return self._last[1]
        N = p * 1e5 / (kB * T)
        env = {"EN": EN, "E": EN * 1e-21 * N, "N": N, "p": p, "T": T, "kB": kB, "qe": Q}
        env.update(self._values)
        for comp, expr in self._gas.items():
            env[self._fraction_names[comp]] = expr(env)
        for name, expr in self._definitions:
            env[name] = expr(env)
        if self._alpha is not None:
            env["alpha"] = self._alpha(env)
        if self._eta is not None:
            env["eta"] = self._eta(env)
        self._last = (key, env)
        return env

    def _density(self, name, env):
        if name == TOTAL_DENSITY:
            return env["N"]
        return env[self._fraction_names[name]] * env["N"]

    def get_R(self, EN, p, T):
        env = self._env(EN, p, T)
        R = np.zeros((len(self.SPECIES), len(self.SPECIES)))
        for i, rx in enumerate(self._reactions):
            K = rx["rate"](env) * self._multipliers.get(i, 1.0)
            for name, count in rx["neutral"].items():
                K *= self._density(name, env) ** count
            R[:, rx["driver"]] += rx["delta"] * K
        return R

    def get_V(self, EN, p, T):
        env = self._env(EN, p, T)
        muN = np.array([m(env) for m in self._mobility])
        return np.diag(self._sign * muN * EN * 1e-21)

    def _selector(self, rows):
        P = np.zeros((len(rows), len(self.SPECIES)))
        for r, i in enumerate(rows):
            P[r, i] = 1.0
        return P

    def get_Pi_e(self):
        return self._selector([self.ELECTRON_INDEX])

    def get_Pi_plus(self):
        return self._selector(self._positive)

    def get_Pi_minus(self):
        return self._selector(self._negative)

    def get_gamma_plus(self, EN, p, T):
        env = self._env(EN, p, T)
        return np.array([g(env) for g in self._gamma_plus])

    def get_gamma_plus_with(
        self, EN, p, T, gamma0=None, gamma1=None, eref=None, beta=None
    ):
        """As :meth:`get_gamma_plus`; per-polarity values come from :meth:`resolve`."""
        return self.get_gamma_plus(EN, p, T)

    def get_B(self, EN, p, T):
        return np.zeros((len(self.SPECIES), 0))

    def get_C(self, EN, p, T):
        return np.zeros((0, len(self.SPECIES)))

    def get_kappa(self, p, T):
        return np.zeros(0)

    def get_gamma_Psi(self, EN, p, T):
        return np.zeros(0)

    # --- polarity ------------------------------------------------------------

    @property
    def has_polarity_overrides(self):
        """True when the configuration sets parameters for one polarity only."""
        return bool(self._polarity)

    def resolve(self, positive):
        """This mechanism with the parameters of one polarity applied."""
        side = "positive" if positive else "negative"
        if side not in self._polarity:
            return self
        config = dict(self._values)
        config.update(self._polarity[side])
        config["label"] = self.label
        resolved = DeclarativeMechanism(self._spec, config, self._source)
        resolved._multipliers = self._multipliers
        return resolved
