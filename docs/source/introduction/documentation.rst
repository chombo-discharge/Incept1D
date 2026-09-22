.. _Chap:Documentation:

Using this documentation
========================

This is the user documentation for ``Incept1D``.  It is organised in five
parts:

* **Introduction** — what you need, where to get the code, and how to run a
  first calculation.
* **Theory** — the one-dimensional drift-reaction model, the two-stream
  photoionization approximation, secondary electron emission, and the
  determinant inception criterion.  Every equation is cross-referenced to
  the function that implements it.
* **Numerics** — how the propagator is built, how the determinant is
  evaluated robustly, how roots are found and tracked, and what a
  calculation costs.
* **Python modules** — what each module and command does, its command-line interface,
  and how to modify or write a reaction mechanism.
* **Examples** — worked reproductions of the comparisons against the IEC
  60052 sphere-gap standard and the Dakin *et al.* (ELECTRA) breakdown
  compilation.
* **Maintenance** — the test suite, continuous integration, and how to
  contribute.

The documentation is built with `Sphinx <https://www.sphinx-doc.org>`_ from
the ``docs/source`` directory of the repository, and the API reference is
generated directly from the NumPy-style docstrings in the code.  The
version at `chombo-discharge.github.io/Incept1D
<https://chombo-discharge.github.io/Incept1D/>`_ is rebuilt by continuous
integration from every commit to ``main``; see :ref:`Chap:Infrastructure`
for how to build it locally.

Notation
--------

Throughout the theory pages, bold upright symbols (:math:`\bm{R}`,
:math:`\bm{V}`, :math:`\bm{Q}`) are matrices and arrow symbols
(:math:`\vec{n}`, :math:`\vec{\theta}`) are column vectors.  Reduced fields :math:`E/N` are given in Townsend
(:math:`1\,\mathrm{Td} = 10^{-21}\,\mathrm{V\,m^2}`), pressures in bar, gap
distances in mm, and the product :math:`pd` in bar·mm.  Inside the code all
quantities are SI except where a function docstring explicitly says
otherwise (the command-line interfaces accept bar, mm, kV and Td).
