.. _Chap:TestSuite:

Test suite
==========

.. note::

   The test suite is not written yet.  This page records the intended
   structure so that tests and documentation are added together.

Planned layout
--------------

Tests will live in ``tests/`` at the repository root and run with
``pytest`` from the repository root:

.. code-block:: bash

   python3 -m pytest tests/

Three layers are planned:

**Unit tests** of the building blocks that have an analytic answer:

* :mod:`incept1d.reactions` — parsing rules, stoichiometry, multiplier
  normalisation, rejection of reactions with two tracked reactants;
* :mod:`incept1d.fields` — normalisation :math:`\int_0^1 f\,d\xi = 1`
  for every geometry, symmetry flags, field-line file layouts and units;
* the propagators — ``magnus2`` reduces to the midpoint rule for a
  constant :math:`\bm{\mathcal{A}}`; composed propagators reproduce a
  single ``expm`` for a uniform field; adaptive refinement respects the
  ``N_max`` budget.

**Verification tests** against closed-form results
(:ref:`Chap:InceptionCriterion`):

* a three-species toy mechanism (ionization, attachment, detachment) must
  reproduce :eq:`eq_generalized_paschen` to solver tolerance, and
  :eq:`eq_standard_paschen` with detachment switched off;
* the leading eigenvalue of the toy mechanism must equal
  :math:`\lambda_+`.

**Example (regression) tests** that run the shipped mechanism on a few
:math:`pd` points of the :ref:`Chap:Examples:IEC60052` and
:ref:`Chap:Examples:Electra` cases with ``--no-plot`` and compare against
stored reference values with a relative tolerance, so that a change in the
chemistry or the numerics that shifts the inception curve is caught.
