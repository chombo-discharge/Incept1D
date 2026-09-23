.. _Chap:Contributing:

How to contribute
=================

.. contents:: On this page
   :local:
   :depth: 1

Contributions — bug reports, new mechanisms, numerical improvements,
documentation fixes — are welcome through the GitHub repository |repo|.

Setting up
----------

.. code-block:: bash

   git clone git@github.com:chombo-discharge/Incept1D.git
   cd Incept1D
   python3 -m venv .venv && source .venv/bin/activate
   pip install numpy scipy matplotlib pre-commit pytest -r docs/requirements.txt
   pre-commit install

``pre-commit install`` registers the hooks that run on every commit:
``black`` (line length 88), ``flake8``, trailing-whitespace and
end-of-file fixes, JSON/YAML validation, a large-file guard (2 MB), and a
Sphinx dummy build.  Run them on demand with

.. code-block:: bash

   pre-commit run --all-files

Raw reference data under ``mechanisms/air/*.txt`` and ``mechanisms/air/*.dat`` (BOLSIG+ output,
LXCat tables, Zheleznyak fit data) is excluded from the whitespace hooks
and must be left byte-for-byte as provided by its source.

Coding conventions
------------------

* Python ≥ 3.10, NumPy/SciPy/matplotlib only; no further runtime
  dependencies.
* Flat layout: modules at the repository root import each other directly
  (``from Constants import kB``).  Do not introduce package-relative
  imports or an ``src/`` layout; mechanism files resolve the root with
  ``sys.path.insert(0, ...)`` relative to their own location.
* Physical constants come from :mod:`incept1d.constants`; never hard-code them.
* Public functions and classes carry NumPy-style docstrings (they are the
  API reference).  Private helpers are prefixed with ``_``.
* Units: SI inside the code; bar / mm / kV / Td on the command line and in
  output files, always labelled.
* ``black`` formats everything; where a file predates ``black`` expect a
  formatting-only diff the first time it is touched.

Physics-affecting changes
-------------------------

If you change a rate coefficient, a boundary condition, a transport
coefficient, or the assembly of :math:`\bm{\mathcal{A}}` or
:math:`\bm{Q}`:

* Cite the equation of :ref:`Chap:TheoryOverview` (or the literature
  source of the new data) in the docstring and the commit message;
* Update the corresponding theory page if the model itself changed;
* Check the closed-form limit with ``mechanisms/air/paschen.json``
  (:eq:`eq_standard_paschen`) and re-run one of the examples;
* Extend the test suite (:ref:`Chap:TestSuite`) where the change has a
  checkable consequence.

New mechanisms
--------------

A new gas goes in its own directory under ``mechanisms/``
(``mechanisms/<gas>/<gas>_<scheme>.py`` + ``mechanisms/<gas>/config.py`` +
data files), following :ref:`Chap:NewMechanisms`.
Please include the swarm-data source, the literature reference for every
rate, and a configuration file that reduces the scheme to the textbook
limit for testing.

Licensing and REUSE
-------------------

The project is licensed under **GPL-3.0-or-later** and follows the
`REUSE <https://reuse.software>`_ specification: *every* file must declare its
copyright holder and licence.  ``reuse lint`` runs both as a pre-commit hook
and in continuous integration, so a file that declares nothing fails the build.

For new source files nothing has to be done by hand beyond the header:

.. code-block:: python

   # SPDX-FileCopyrightText: 2026 SINTEF Energy Research
   #
   # SPDX-License-Identifier: GPL-3.0-or-later

Files where a comment header would be rendered by the consumer (``.rst``,
``.json``, ``.css``) are declared in bulk in ``REUSE.toml`` instead; the
existing glob patterns already cover new files of those types.

Third-party data
................

Some data in this repository is not ours.  Redistributed swarm data, for
instance, belongs to the databases it was retrieved from, is *not* under the
project licence, and carries its real rights holder and a ``LicenseRef-``
licence in ``REUSE.toml``, with the terms spelled out in ``LICENSES/``.

.. warning::

   **Do not add content that this project has no right to redistribute.**
   Published standards, journal tables, figures and datasets are copyrighted
   by their publishers, and a small extract used for validation is still
   redistribution.  If you cannot point to a licence or a permission that
   allows it, it does not go in the repository — not in ``examples/``, not in
   a figure, not in a docstring.

   This is not a formality.  A comparison against data you may not
   redistribute is still perfectly reproducible: commit the *calculation*,
   document the expected input file and its column layout, and let a reader
   who holds the data supply it themselves.  Two of the worked examples are
   built exactly that way, and the figures they would produce are opt-in for
   the same reason (:ref:`Chap:Infrastructure`).

When you add data that *is* redistributable:

* Redistribute it **verbatim**, with its original header intact — the header
  carries the citation its authors ask for, and the whitespace hooks skip
  ``mechanisms/**`` and ``examples/**`` data files on purpose.
* Add a ``[[annotations]]`` entry naming the actual rights holder, never
  ``SINTEF Energy Research`` by default.
* If its terms are not an existing SPDX licence, add a
  ``LICENSES/LicenseRef-<name>.txt`` describing them, and reference it.
* Record where it came from and when, so the claim can be checked later.

Mislabelling third-party data as project-owned is worse than leaving it
undeclared, because ``reuse lint`` will then happily pass.

Documentation
-------------

Documentation lives in ``docs/source`` and is built as described in
:ref:`Chap:Infrastructure`.  New pages are added to the appropriate
``toctree`` in ``index.rst``.  A pull request that changes a public
function should update its docstring; a pull request that changes the
command-line interface should update the corresponding page under
*Python modules*.

Submitting changes
------------------

1. Create a branch from ``main``.
2. Make the change with tests and documentation.
3. Ensure ``pre-commit run --all-files`` and ``python3 -m pytest`` pass.
4. Open a pull request describing *what* changed and *why*; for physics
   changes include a before/after inception curve or equivalent evidence.
