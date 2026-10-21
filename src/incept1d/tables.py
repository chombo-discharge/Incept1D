# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Tabulated swarm data for declarative mechanism files.

Two formats are read, each into columns tabulated against E/N [Td]:

``bolsig+``
    BOLSIG+ output.  Section 1 holds the transport coefficients, one column
    per label (``A1`` mean energy [eV], ``A2`` mobility * N [1/(m V s)],
    ``A6`` diffusion * N [1/(m s)], ...); section 2 the rate coefficient of
    each collision process [m^3/s], named in its legend by target and type
    (``C25  N2  Ionization  15.60 eV``).  A rate is selected by target and
    type, ``"N2 Ionization"``, summing every process of that kind.  The
    mixture (``Mole fraction N2  0.7900``) and gas temperature (``Gas
    temperature (K)`` or ``Tgas = 300 K``) are read when the file records
    them.
``lxcat``
    A two-column LXCat table (E/N [Td], value) between the first pair of
    ``----`` lines, such as the reduced ion mobilities K0 [cm^2/(V s)] of the
    Viehland database; column 1 is the value.

Each file is parsed once per process.
"""

import os
import re

import numpy as np

_CACHE = {}


class TableError(ValueError):
    """A table that cannot be read, or a column it does not have."""


class Table:
    """
    Columns tabulated against E/N.

    Attributes
    ----------
    path : str
    fractions : dict
        {gas component: mole fraction} recorded in the file, or empty.
    T_gas : float or None
        Gas temperature [K] recorded in the file.
    """

    def __init__(self, path, fmt):
        self.path = path
        self.format = fmt
        self.fractions = {}
        self.T_gas = None
        self._columns = {}  # label -> (EN, values)
        self._rates = {}  # (target, type) -> list of (EN, values)

    def column(self, label):
        """(E/N, values) of a column, by its label."""
        key = str(label)
        if key not in self._columns:
            raise TableError(
                f"{os.path.basename(self.path)} has no column {label!r}; "
                f"it has {sorted(self._columns)}"
            )
        return self._columns[key]

    def rate(self, process):
        """(E/N, summed rate coefficient) of every process 'target type'."""
        parts = str(process).split()
        if self.format != "bolsig+" or len(parts) != 2:
            raise TableError(
                f"a rate is selected as 'target type', e.g. 'N2 Ionization', "
                f"from a bolsig+ table; got {process!r}"
            )
        found = self._rates.get(tuple(parts))
        if not found:
            raise TableError(
                f"{os.path.basename(self.path)} lists no {process!r} process; "
                f"it lists {sorted(' '.join(k) for k in self._rates)}"
            )
        EN = found[0][0]
        total = np.zeros(len(EN))
        for _, values in found:
            total += values
        return EN, total


def read_table(path, fmt):
    """Parse *path* in format *fmt* (``"bolsig+"`` or ``"lxcat"``), once."""
    key = (os.path.realpath(path), fmt)
    if key not in _CACHE:
        if not os.path.isfile(path):
            raise TableError(f"table file not found: {path}")
        readers = {"bolsig+": _read_bolsig, "lxcat": _read_lxcat}
        if fmt not in readers:
            raise TableError(f"unknown table format {fmt!r}; known: {sorted(readers)}")
        table = Table(path, fmt)
        readers[fmt](path, table)
        _CACHE[key] = table
    return _CACHE[key]


def _rows(path, lines, start):
    """Numeric rows from *start* up to the first blank or non-numeric line."""
    rows = []
    for number, ln in enumerate(lines[start:], start=start + 1):
        vals = ln.split()
        if not vals or not vals[0].lstrip("-").isdigit():
            break
        try:
            rows.append([float(v) for v in vals])
        except ValueError:
            raise TableError(f"{path}: line {number} is not numeric: {ln!r}") from None
    if len({len(r) for r in rows}) > 1:
        raise TableError(f"{path}: rows of unequal length after line {start}")
    return np.array(rows)


def _check_width(path, section, data, columns):
    """Raise TableError unless *data* has rows with at least *columns* columns."""
    if data.ndim != 2 or len(data) == 0:
        raise TableError(f"{path}: the {section} section has no data rows")
    if data.shape[1] < columns:
        raise TableError(
            f"{path}: the {section} section has {data.shape[1]} columns, but its "
            f"header names {columns}"
        )


def _read_bolsig(path, table):
    with open(path, "rb") as fh:
        text = fh.read().decode("ascii", errors="replace")
    lines = [ln.strip() for ln in text.splitlines()]
    table.fractions = {
        m.group(1): float(m.group(2))
        for m in re.finditer(r"Mole fraction\s+(\S+)\s+([-+0-9.Ee]+)", text)
    }
    m = re.search(r"Gas temperature \(K\)\s+([-+0-9.Ee]+)", text) or re.search(
        r"Tgas\s*=\s*([-+0-9.Ee]+)\s*K", text
    )
    table.T_gas = float(m.group(1)) if m else None

    # Section 1: "R#  E/N (Td)  A1  A2  A6 ...", data columns R#, E/N, then the labels.
    head = next(
        (
            i
            for i, ln in enumerate(lines)
            if ln.startswith("R#")
            and "E/N (Td)" in ln
            and "A1" in ln
            and "Energy" not in ln
        ),
        None,
    )
    if head is None:
        raise TableError(f"{path}: no transport section ('R#  E/N (Td)  A1 ...')")
    labels = lines[head].split("E/N (Td)", 1)[1].split()
    data = _rows(path, lines, head + 1)
    _check_width(path, "transport", data, 2 + len(labels))
    for k, label in enumerate(labels):
        table._columns[label] = (data[:, 1], data[:, 2 + k])

    # Section 2: legend "Cx  target  type ...", then "R#  E/N (Td)  Energy (eV)  C1 ...".
    try:
        start = next(
            i
            for i, ln in enumerate(lines)
            if "Rate coefficients (m3/s)" in ln and "Inverse" not in ln
        )
    except StopIteration:
        return
    legend = []
    for i in range(start + 1, len(lines)):
        ln = lines[i]
        if ln.startswith("R#") and "E/N" in ln:
            data = _rows(path, lines, i + 1)
            _check_width(path, "rate coefficient", data, 3)
            break
        parts = ln.split()
        if parts and parts[0].startswith("C") and parts[0][1:].isdigit():
            legend.append((int(parts[0][1:]), tuple(parts[1:3])))
    else:
        return
    for cx, kind in legend:
        if len(kind) == 2 and cx + 2 < data.shape[1]:
            table._rates.setdefault(kind, []).append((data[:, 1], data[:, cx + 2]))


def _read_lxcat(path, table):
    rows, inside = [], False
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("----"):
                if inside:
                    break
                inside = True
                continue
            if inside:
                parts = stripped.split()
                if len(parts) == 2:
                    rows.append([float(parts[0]), float(parts[1])])
    if not rows:
        raise TableError(f"{path}: no data between '----' lines")
    data = np.array(rows)
    table._columns["1"] = (data[:, 0], data[:, 1])
