# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Declarative mechanism files: a JSON description compiled into the solver
interface.

A mechanism file names the species, their charges and mobilities, the
reactions with their rate coefficients, the Townsend coefficients, the
photoionization model and the electrode yields, with every coefficient a
formula (:class:`incept1d.expressions.Expression`) over named variables.
:func:`load_json_mechanism` checks the file strictly — an unknown key, an
unknown name in a formula or a reaction that does not balance is an error,
not a silent default — and returns a :class:`DeclarativeMechanism`, which
has the interface the solver uses (see
:data:`incept1d.mechanism.REQUIRED_ATTRS`).  What the Python interface asks
a mechanism to build by hand follows from the file: the row selectors and
the drift directions from the charges, the reaction matrix from the
reaction strings, the photon matrices from the photoionization block.

A configuration (``read_json_configs``) sets the file's parameters by name,
scales reactions by their reaction string (``reaction_multipliers``) and
overrides parameters for one polarity (``positive``, ``negative``).  A
parameter whose value is a string names a file; a relative path in the
mechanism file is relative to it, one in a configuration to the
configuration file, else to the mechanism file.

Formulas are evaluated with these variables:

    EN   reduced field [Td]          E   field [V/m]
    N    neutral number density [m^-3]
    p    pressure [bar]              T   gas temperature [K]
    kB   Boltzmann constant [J/K]    qe  elementary charge [C]
    pi
    x_<X>  mole fraction of gas component X
    every numeric parameter, every table column and every definition
    alpha, eta (when the file gives them), after the definitions
"""

import json
import os
import re
import warnings

import numpy as np

from incept1d.constants import Q, kB
from incept1d.expressions import FUNCTIONS, Expression
from incept1d.reactions import _normalize, parse_reaction
from incept1d.tables import read_table

#: The ``format`` a mechanism file must declare.
FORMAT = "incept1d-mechanism/1"

_TOP_KEYS = {
    "format",
    "description",
    "parameters",
    "gas",
    "tables",
    "untracked",
    "definitions",
    "species",
    "alpha",
    "eta",
    "reactions",
    "photoionization",
    "emission",
}
_SPECIES_KEYS = {"name", "charge", "mobility"}
_REACTION_KEYS = {"reaction", "rate"}
_TABLE_KEYS = {"format", "file", "columns", "rates"}
_PHOTO_KEYS = {"absorption", "source", "products", "cathode yield", "direction"}
_ABSORPTION_KEYS = {"model", "groups", "absorber"}
_SOURCE_KEYS = {"reactions", "efficiency"}
_EMISSION_KEYS = {"ions"}
_CONFIG_KEYS = {"label", "reaction_multipliers", "positive", "negative"}
_BUILTINS = ("EN", "E", "N", "p", "T", "kB", "qe", "pi")
_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_DIRECTIONS = ("both", "forward", "backward")

#: The tracked species that is the electron.
ELECTRON = "e"
#: The name that stands for the total neutral density in a reaction.
TOTAL_DENSITY = "M"
#: Largest difference between the gas mixture and a table's recorded one.
MIXTURE_TOL = 1e-3
#: Gas temperature difference [K] from a table's above which a warning is given.
TEMPERATURE_TOL = 20.0

_WARNED = set()


class MechanismFileError(ValueError):
    """A mechanism file or configuration that does not describe a mechanism."""


def is_mechanism_file(path):
    """True if *path* is a JSON file declaring the mechanism :data:`FORMAT`."""
    if not str(path).endswith(".json"):
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
    return DeclarativeMechanism(spec, config_dict, source=os.path.abspath(path))


def _keys(obj, allowed, where):
    """Check that *obj* is a dict whose keys are *allowed* (or start with _)."""
    if not isinstance(obj, dict):
        raise MechanismFileError(f"{where}: expected an object, got {obj!r}")
    unknown = sorted(k for k in obj if k not in allowed and not k.startswith("_"))
    if unknown:
        raise MechanismFileError(
            f"{where}: unknown key(s) {unknown}; expected some of {sorted(allowed)}"
        )


def _require(obj, keys, where):
    for k in keys:
        if k not in obj:
            raise MechanismFileError(f"{where}: missing {k!r}")


def _items(obj):
    """The entries of a JSON object that are not comments."""
    return [(k, v) for k, v in obj.items() if not k.startswith("_")]


def _warn_once(key, message):
    if key not in _WARNED:
        _WARNED.add(key)
        warnings.warn(message, stacklevel=4)


class DeclarativeMechanism:
    """
    A mechanism compiled from a declarative description.

    Parameters
    ----------
    spec : dict
        The parsed mechanism file.
    config : dict or None
        One configuration: parameter values, ``label``,
        ``reaction_multipliers``, ``positive`` / ``negative``, and
        ``_cfg_dir``, the directory relative file names in it are taken from.
    source : str
        The mechanism file, for error messages and relative file names.
    """

    def __init__(self, spec, config=None, source="<mechanism>"):
        self._spec = spec
        self._source = source
        self._dir = os.path.dirname(source) if os.path.isfile(source) else os.getcwd()
        _keys(spec, _TOP_KEYS, source)
        if spec.get("format") != FORMAT:
            raise MechanismFileError(
                f"{source}: 'format' must be {FORMAT!r}, got {spec.get('format')!r}"
            )
        self._parameters = self._read_parameters(spec.get("parameters", {}))
        config = dict(config or {})
        _keys(config, _CONFIG_KEYS | set(self._parameters), f"{source}: configuration")
        self._cfg_dir = config.get("_cfg_dir", self._dir)
        self.label = str(config.get("label", "Baseline"))
        self._configured = {k for k in config if k in self._parameters}
        self._values = dict(self._parameters)
        self._values.update(self._overrides(config, f"{source}: configuration"))
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
        for name, value in _items(params):
            self._check_name(name, "parameter")
            where = f"{self._source}: parameter {name!r}"
            if isinstance(value, str):
                out[name] = self._path(value, self._dir)
            elif isinstance(value, bool) or not isinstance(value, (int, float)):
                raise MechanismFileError(
                    f"{where}: expected a number or a file name, got {value!r}"
                )
            else:
                out[name] = float(value)
        return out

    def _overrides(self, values, where):
        """The parameter values in *values*, of the parameters' own types."""
        out = {}
        for k, v in values.items():
            if k not in self._parameters:
                continue
            if isinstance(self._parameters[k], str):
                if not isinstance(v, str):
                    raise MechanismFileError(f"{where} {k!r}: expected a file name")
                candidate = self._path(v, self._cfg_dir)
                out[k] = (
                    candidate if os.path.isfile(candidate) else self._path(v, self._dir)
                )
            elif isinstance(v, bool) or not isinstance(v, (int, float)):
                raise MechanismFileError(f"{where} {k!r}: expected a number, got {v!r}")
            else:
                out[k] = float(v)
        return out

    @staticmethod
    def _path(name, base):
        return (
            name if os.path.isabs(name) else os.path.normpath(os.path.join(base, name))
        )

    def _check_name(self, name, what):
        if not _NAME.match(name) or name in _BUILTINS or name in FUNCTIONS:
            raise MechanismFileError(
                f"{self._source}: {what} name {name!r} must start with a letter, "
                f"contain only letters, digits and _, and not be one of "
                f"{sorted(set(_BUILTINS) | set(FUNCTIONS))}"
            )

    def _declare(self, name, names, what):
        self._check_name(name, what)
        if name in names or name in self._parameters:
            raise MechanismFileError(
                f"{self._source}: {what} {name!r} is already defined"
            )
        names.append(name)

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
            elif isinstance(value, bool) or not isinstance(value, (int, float)):
                raise MechanismFileError(
                    f"{self._source}: multiplier {key!r}: not a number"
                )
            else:
                out[index] = float(value)
        if unknown:
            raise MechanismFileError(
                f"Configuration '{self.label}' sets multipliers for reactions that "
                f"'{os.path.basename(self._source)}' does not define: {unknown}"
            )
        return out

    # --- compiling -----------------------------------------------------------

    def _compile(self):
        spec, src = self._spec, self._source
        numeric = [k for k, v in self._parameters.items() if not isinstance(v, str)]
        names = list(_BUILTINS) + numeric
        self._constants = {"kB": kB, "qe": Q, "pi": np.pi}
        self._constants.update({k: v for k, v in self._values.items() if k in numeric})

        self._compile_gas(spec.get("gas", {}), numeric, names)
        self._compile_tables(spec.get("tables", {}), names)

        self._definitions = []
        defs = spec.get("definitions") or {}
        _keys(defs, {k for k in defs}, f"{src}: definitions")
        for name, text in _items(defs):
            expr = Expression(text, names, f"{src}: definition {name!r}")
            self._declare(name, names, "definition")
            self._definitions.append((name, expr))

        self._alpha = self._eta = None
        if "alpha" in spec:
            self._alpha = Expression(spec["alpha"], names, f"{src}: alpha")
        if "eta" in spec:
            self._eta = Expression(spec["eta"], names, f"{src}: eta")
        names += [n for n, e in (("alpha", self._alpha), ("eta", self._eta)) if e]
        if self._alpha is not None:
            self.alpha = lambda EN, p, T: self._env(EN, p, T)["alpha"]
        if self._eta is not None:
            self.eta = lambda EN, p, T: self._env(EN, p, T)["eta"]

        self._compile_species(spec.get("species"), names)

        untracked = spec.get("untracked", [])
        if not isinstance(untracked, list) or not all(
            isinstance(u, str) for u in untracked
        ):
            raise MechanismFileError(f"{src}: 'untracked' must be a list of names")
        reactions = spec.get("reactions") or []
        if not isinstance(reactions, list):
            raise MechanismFileError(f"{src}: 'reactions' must be a list")
        self._reactions = []
        seen = set()
        for i, rx in enumerate(reactions):
            where = f"{src}: reactions[{i}]"
            _keys(rx, _REACTION_KEYS, where)
            _require(rx, _REACTION_KEYS, where)
            compiled = self._compile_reaction(rx, set(untracked), names, where)
            key = _normalize(compiled["reaction"])
            if key in seen:
                raise MechanismFileError(
                    f"{where}: {compiled['reaction']!r} appears twice; the reaction "
                    f"string is the reaction's name"
                )
            seen.add(key)
            self._reactions.append(compiled)

        self._compile_photoionization(spec.get("photoionization"), names)
        self._compile_emission(spec.get("emission", {}), names)

    def _compile_gas(self, gas, numeric, names):
        src = self._source
        if not isinstance(gas, dict):
            raise MechanismFileError(f"{src}: 'gas' must be an object")
        self._gas = {}
        for comp, frac in _items(gas):
            self._gas[comp] = Expression(frac, numeric, f"{src}: gas {comp!r}")
        self._fractions = {c: e(self._constants) for c, e in self._gas.items()}
        if self._fractions:
            bad = {c: x for c, x in self._fractions.items() if not 0.0 <= x <= 1.0}
            if bad or abs(sum(self._fractions.values()) - 1.0) > MIXTURE_TOL:
                raise MechanismFileError(
                    f"{src}: gas mole fractions {self._fractions} must lie in [0, 1] "
                    f"and sum to 1"
                )
        self._fraction_names = {}
        for comp in self._gas:
            self._declare(f"x_{comp}", names, f"mole fraction of {comp!r}")
            self._fraction_names[comp] = f"x_{comp}"
        self._constants.update(
            {self._fraction_names[c]: x for c, x in self._fractions.items()}
        )
        # Parameters that set the mixture, to tell a configured one from defaults.
        self._gas_parameters = set().union(*(e.names for e in self._gas.values()))

    def _compile_tables(self, tables, names):
        src = self._source
        if not isinstance(tables, dict):
            raise MechanismFileError(f"{src}: 'tables' must be an object")
        self._columns = []  # (variable, EN, values)
        self._checked_tables = []
        for tname, tspec in _items(tables):
            where = f"{src}: tables[{tname!r}]"
            _keys(tspec, _TABLE_KEYS, where)
            _require(tspec, ("format", "file"), where)
            file = tspec["file"]
            if isinstance(self._parameters.get(file), str):
                path = self._values[file]
            elif isinstance(file, str):
                path = self._path(file, self._dir)
            else:
                raise MechanismFileError(
                    f"{where}: 'file' must be a file or parameter name"
                )
            try:
                table = read_table(path, tspec["format"])
                for kind in ("columns", "rates"):
                    entries = tspec.get(kind, {})
                    _keys(entries, {k for k in entries}, f"{where}.{kind}")
                    for var, label in _items(entries):
                        x, y = (
                            table.column(label)
                            if kind == "columns"
                            else table.rate(label)
                        )
                        self._declare(var, names, "table column")
                        self._columns.append((var, x, y))
            except ValueError as exc:
                if isinstance(exc, MechanismFileError):
                    raise
                raise MechanismFileError(f"{where}: {exc}") from None
            if table.format == "bolsig+":
                self._check_mixture(table, where)
                self._checked_tables.append(table)

    def _check_mixture(self, table, where):
        """The neutral densities must be those of the mixture the table is for."""
        name = os.path.basename(table.path)
        if table.fractions:
            for comp in set(self._fractions) | set(table.fractions):
                ours = self._fractions.get(comp, 0.0)
                theirs = table.fractions.get(comp, 0.0)
                if abs(ours - theirs) > MIXTURE_TOL:
                    raise MechanismFileError(
                        f"{where}: the gas mixture {self._fractions} differs from the "
                        f"one {name} was computed for, {table.fractions}.  Use "
                        f"transport data for the mixture, or set the mixture."
                    )
        elif self._fractions and not (self._gas_parameters & self._configured):
            _warn_once(
                ("mixture", table.path, tuple(sorted(self._fractions.items()))),
                f"{name} does not record the gas mixture of its transport data; "
                f"using {self._fractions}.  Set the mixture the table was computed "
                f"for in the configuration.",
            )

    def _compile_species(self, species, names):
        src = self._source
        if not isinstance(species, list) or not species:
            raise MechanismFileError(f"{src}: 'species' must be a non-empty list")
        self.SPECIES, self._charge, self._mobility = [], [], []
        for i, sp in enumerate(species):
            where = f"{src}: species[{i}]"
            _keys(sp, _SPECIES_KEYS, where)
            _require(sp, _SPECIES_KEYS, where)
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
            if not (
                name in self.SPECIES
                or name in self._gas
                or name == TOTAL_DENSITY
                or name in untracked
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

    def _compile_photoionization(self, photo, names):
        src = self._source
        self._n_gamma = 0
        self._photo = None
        self.PHOTON_DIRECTIONS = None
        if photo is None:
            return
        where = f"{src}: photoionization"
        _keys(photo, _PHOTO_KEYS, where)
        _require(photo, ("absorption", "source", "products", "cathode yield"), where)

        absorption = photo["absorption"]
        _keys(absorption, _ABSORPTION_KEYS, f"{where}.absorption")
        _require(absorption, _ABSORPTION_KEYS, f"{where}.absorption")
        if absorption["model"] != "zheleznyak":
            raise MechanismFileError(
                f"{where}.absorption: unknown model {absorption['model']!r}; "
                f"known: ['zheleznyak']"
            )
        numeric = [k for k, v in self._parameters.items() if not isinstance(v, str)]
        groups = Expression(absorption["groups"], numeric, f"{where}.absorption.groups")
        n = groups(self._constants)
        if n != int(n) or n < 1:
            raise MechanismFileError(
                f"{where}.absorption: 'groups' = {n:g} is not a count"
            )
        absorber = absorption["absorber"]
        if absorber not in self._gas:
            raise MechanismFileError(
                f"{where}.absorption: absorber {absorber!r} is not a gas component"
            )
        from incept1d.zheleznyak import twostream_groups

        kappa, g = twostream_groups(int(n))

        source = photo["source"]
        _keys(source, _SOURCE_KEYS, f"{where}.source")
        _require(source, _SOURCE_KEYS, f"{where}.source")
        known = {_normalize(r["reaction"]): r for r in self._reactions}
        if not isinstance(source["reactions"], list) or not source["reactions"]:
            raise MechanismFileError(
                f"{where}.source: 'reactions' must be a non-empty list"
            )
        emitting = []
        for text in source["reactions"]:
            rx = known.get(_normalize(str(text)))
            if rx is None:
                raise MechanismFileError(f"{where}.source: no reaction {text!r}")
            emitting.append(rx)

        products = photo["products"]
        _keys(products, set(self.SPECIES), f"{where}.products")
        direction = photo.get("direction", "both")
        if direction not in _DIRECTIONS:
            raise MechanismFileError(f"{where}: direction must be one of {_DIRECTIONS}")
        self._n_gamma = int(n)
        if "direction" in photo:
            self.PHOTON_DIRECTIONS = [direction] * self._n_gamma
        self._photo = {
            "kappa": np.array(kappa) * self._fractions[absorber] * 1e5,  # per bar of p
            "g": np.array(g),
            "emitting": emitting,
            "efficiency": Expression(
                source["efficiency"], names, f"{where}.source.efficiency"
            ),
            "products": [
                (
                    self.SPECIES.index(s),
                    Expression(e, names, f"{where}.products[{s!r}]"),
                )
                for s, e in _items(products)
            ],
            "yield": Expression(
                photo["cathode yield"], names, f"{where}.cathode yield"
            ),
        }

    def _compile_emission(self, emission, names):
        src = self._source
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

    # --- evaluating ----------------------------------------------------------

    def _env(self, EN, p, T):
        """All variables at (EN, p, T); the last one is kept."""
        key = (EN, p, T)
        if self._last is not None and self._last[0] == key:
            return self._last[1]
        for table in self._checked_tables:
            if table.T_gas is not None and abs(T - table.T_gas) > TEMPERATURE_TOL:
                _warn_once(
                    ("temperature", table.path, float(T)),
                    f"Gas temperature T = {T:g} K differs from the {table.T_gas:g} K "
                    f"of the transport data in {os.path.basename(table.path)}.",
                )
        N = p * 1e5 / (kB * T)
        env = dict(self._constants)
        env.update({"EN": EN, "E": EN * 1e-21 * N, "N": N, "p": p, "T": T})
        for var, x, y in self._columns:
            env[var] = float(np.interp(EN, x, y))
        for name, expr in self._definitions:
            env[name] = expr(env)
        if self._alpha is not None:
            env["alpha"] = self._alpha(env)
        if self._eta is not None:
            env["eta"] = self._eta(env)
        self._last = (key, env)
        return env

    def _rate(self, rx, env):
        """K [1/s] of a reaction, before multipliers."""
        K = rx["rate"](env)
        for name, count in rx["neutral"].items():
            density = (
                env["N"]
                if name == TOTAL_DENSITY
                else (env[self._fraction_names[name]] * env["N"])
            )
            K *= density**count
        return K

    def get_R(self, EN, p, T):
        env = self._env(EN, p, T)
        R = np.zeros((len(self.SPECIES), len(self.SPECIES)))
        for i, rx in enumerate(self._reactions):
            K = self._rate(rx, env) * self._multipliers.get(i, 1.0)
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

    # --- photons -------------------------------------------------------------

    def get_kappa(self, p, T):
        if self._photo is None:
            return np.zeros(0)
        return self._photo["kappa"] * p

    def get_B(self, EN, p, T):
        B = np.zeros((len(self.SPECIES), self._n_gamma))
        if self._photo is not None:
            env = self._env(EN, p, T)
            kappa = self.get_kappa(p, T)
            for i, coeff in self._photo["products"]:
                B[i, :] = coeff(env) * kappa
        return B

    def get_C(self, EN, p, T):
        C = np.zeros((self._n_gamma, len(self.SPECIES)))
        if self._photo is not None:
            env = self._env(EN, p, T)
            emitted = {}
            for rx in self._photo["emitting"]:
                emitted[rx["driver"]] = emitted.get(rx["driver"], 0.0) + self._rate(
                    rx, env
                )
            efficiency = self._photo["efficiency"](env)
            for j, K in emitted.items():
                C[:, j] = efficiency * K * self._photo["g"]
        return C

    def get_gamma_Psi(self, EN, p, T):
        if self._photo is None:
            return np.zeros(0)
        return np.full(self._n_gamma, self._photo["yield"](self._env(EN, p, T)))

    # --- inspection (for exports and diagnostics) ----------------------------

    @property
    def reaction_strings(self):
        """The reactions, as written, in file order."""
        return [rx["reaction"] for rx in self._reactions]

    @property
    def mole_fractions(self):
        """{gas component: mole fraction}."""
        return dict(self._fractions)

    @property
    def files(self):
        """{parameter: resolved path} of the file parameters, as configured."""
        return {k: v for k, v in self._values.items() if isinstance(v, str)}

    def neutral_partners(self, index):
        """{neutral reactant: multiplicity} of reaction *index* ('M' = all)."""
        return dict(self._reactions[index]["neutral"])

    def rate_coefficient(self, index, EN, p, T):
        """Rate coefficient k of reaction *index*, SI, with its multiplier."""
        rx = self._reactions[index]
        return rx["rate"](self._env(EN, p, T)) * self._multipliers.get(index, 1.0)

    def neutral_density_factor(self, index, EN, p, T):
        """The product of the neutral densities reaction *index* is multiplied by."""
        env = self._env(EN, p, T)
        factor = 1.0
        for name, count in self._reactions[index]["neutral"].items():
            n = (
                env["N"]
                if name == TOTAL_DENSITY
                else env[self._fraction_names[name]] * env["N"]
            )
            factor *= n**count
        return factor

    def variable(self, name, EN, p, T):
        """The value of a variable of the file (parameter, column, definition)."""
        env = self._env(EN, p, T)
        if name not in env:
            raise MechanismFileError(f"{self._source} defines no variable {name!r}")
        return env[name]

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
        config = {k: self._values[k] for k in self._configured}
        config.update(self._polarity[side])
        config["label"] = self.label
        resolved = DeclarativeMechanism(self._spec, config, self._source)
        resolved._multipliers = self._multipliers
        return resolved
