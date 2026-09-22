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
* a LaTeX distribution with ``pdflatex`` if you want the PDF version.

Development
-----------

Contributors should also install

* **pre-commit**, which runs ``black``, ``flake8``, and a Sphinx dummy build
  before every commit (see :ref:`Chap:Contributing`);
* **pytest** for the test suite (see :ref:`Chap:TestSuite`).

Input data
----------

The mechanism files shipped in ``mechanisms/air/`` read electron swarm data from
BOLSIG+ output files (``mechanisms/air/lisbon.txt``, ``mechanisms/air/phelps.txt``,
``mechanisms/air/biagi.txt``, ``mechanisms/air/trinity.txt``, ``mechanisms/air/morgan.txt``) and negative-ion
mobility tables from LXCat (``mechanisms/air/o2m_mobility.txt``,
``mechanisms/air/o3m_mobility.txt``).  These are included in the repository; no
external database access is needed to run the examples.  If you build a
mechanism for another gas you will need to generate the corresponding swarm
data yourself with `BOLSIG+ <https://www.bolsig.laplace.univ-tlse.fr>`_ or
an equivalent Boltzmann solver (see :ref:`Chap:NewMechanisms`).
