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

   git clone git@github.com:rmrsk/Incept1D.git
   cd Incept1D
   python3 -m venv .venv && source .venv/bin/activate
   pip install numpy scipy matplotlib pre-commit pytest -r Docs/requirements.txt
   pre-commit install

``pre-commit install`` registers the hooks that run on every commit:
``black`` (line length 88), ``flake8``, trailing-whitespace and
end-of-file fixes, JSON/YAML validation, a large-file guard (2 MB), and a
Sphinx dummy build.  Run them on demand with

.. code-block:: bash

   pre-commit run --all-files

Raw reference data under ``Air/*.txt`` and ``Air/*.dat`` (BOLSIG+ output,
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
* Physical constants come from :mod:`Constants`; never hard-code them.
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

* cite the equation of :ref:`Chap:TheoryOverview` (or the literature
  source of the new data) in the docstring and the commit message;
* update the corresponding theory page if the model itself changed;
* check the closed-form limit with ``Air/Paschen.json``
  (:eq:`eq_standard_paschen`) and re-run one of the examples;
* extend the test suite (:ref:`Chap:TestSuite`) where the change has a
  checkable consequence.

New mechanisms
--------------

A new gas goes in its own directory (``Gas/Gas_<scheme>.py`` +
``Gas/Config.py`` + data files), following :ref:`Chap:NewMechanisms`.
Please include the swarm-data source, the literature reference for every
rate, and a configuration file that reduces the scheme to the textbook
limit for testing.

Documentation
-------------

Documentation lives in ``Docs/source`` and is built as described in
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
   changes include a before/after Paschen curve or equivalent evidence.
