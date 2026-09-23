.. _Chap:Prerequisites:

Prerequisites
=============

``Incept1D`` is pure Python and runs on Linux, macOS, and Windows.  It has
no build step and no compiled components.

Runtime
-------

* **Python** ≥ 3.10.
* **NumPy** ≥ 1.25 (NumPy 2.x is supported).
* **SciPy** ≥ 1.10 — used for ``scipy.linalg.expm`` (matrix exponentials),
  ``scipy.optimize.brentq`` (root finding), ``scipy.optimize.nnls`` (the
  Zheleznyak photoionization fit), ``scipy.optimize.linear_sum_assignment``
  (eigenvalue tracking), and ``scipy.constants``.
* **Matplotlib** ≥ 3.6 — every subcommand produces a figure by
  default.  Headless machines should set ``MPLBACKEND=Agg`` and use the
  ``--no-plot`` flag where available.

Documentation
-------------

Building the documentation additionally requires

* **Sphinx** ≥ 7.0 and **sphinx-rtd-theme** ≥ 2.0 (listed in
  ``docs/requirements.txt``);
* A LaTeX distribution with ``pdflatex`` if you want the PDF version.

Development
-----------

Contributors should also install

* **pre-commit**, which runs ``black``, ``flake8``, and a Sphinx dummy build
  before every commit (see :ref:`Chap:Contributing`);
* **pytest** for the test suite (see :ref:`Chap:TestSuite`).

Input data
----------

The solver itself reads nothing from disk: every gas-dependent quantity comes
from a mechanism file, and what that file reads is its own business.  In
practice a mechanism needs electron swarm data — the ionization, attachment
and transport coefficients as functions of :math:`E/N` — usually produced by a
Boltzmann solver such as `BOLSIG+
<https://www.bolsig.laplace.univ-tlse.fr>`_ from a cross-section set, together
with ion mobilities.

Everything the shipped dry-air mechanism needs is already in the repository,
so no external database access is required to reproduce the examples.  If you
build a mechanism for another gas you will have to supply the equivalent data
yourself; :ref:`Chap:NewMechanisms` describes what is required and how a
mechanism is expected to load it.
