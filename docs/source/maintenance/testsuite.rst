.. _Chap:TestSuite:

Test suite
==========

Layout
------

Tests live in ``tests/`` at the repository root and run with ``pytest``:

.. code-block:: bash

   python3 -m pytest tests/

Install the development extra first (``pip install -e '.[dev]'``).  The
default run takes a few seconds; the tests that load real swarm data are
marked ``slow`` and can be deselected:

.. code-block:: bash

   python3 -m pytest -m "not slow"

Three layers, of decreasing strength:

**Verification against closed-form results** — the core of the suite.
``tests/toy/toy_mechanism.py`` is a three-species mechanism (ionization,
attachment, detachment, cathode secondary emission) whose inception
condition is :eq:`eq_generalized_paschen`.  That algebra is transcribed
independently in ``tests/closed_form.py``, so the two never share code.
The solver must reproduce it, and its documented limits, to about
:math:`10^{-7}` relative:

* :eq:`eq_generalized_paschen` over a range of gap lengths;
* :eq:`eq_standard_paschen` with detachment switched off;
* :math:`\alpha d = \ln(1 + \gamma^{-1})` with attachment switched off too;
* Solutions in the attachment regime :math:`\alpha \le \eta`, where the
  long-gap limit is :math:`\alpha \to \eta/(1+\gamma)`;
* The leading eigenvalue of :math:`\bm{R}\bm{V}^{-1}` equals
  :math:`\lambda_+`.

The shipped dry-air mechanism is held to the same standard:
``mechanisms/air/paschen.json`` reduces it to the textbook limit, and it
must then satisfy :eq:`eq_standard_paschen` using its own tabulated
:math:`\alpha`, :math:`\eta` and :math:`\gamma`.

**Unit tests** of the building blocks:

* :mod:`incept1d.reactions` — parsing rules, stoichiometry, multiplier
  normalisation, rejection of reactions with two tracked reactants;
* :mod:`incept1d.fields` — normalisation :math:`\int_0^1 f\,d\xi = 1`
  for every geometry, symmetry flags, field-line file layouts, units and
  malformed input;
* :mod:`incept1d.mechanism` — the interface contract, the ``config.py``
  protocol and JSON configuration shapes;
* The propagators — ``magnus2`` reduces to the midpoint rule for constant
  :math:`\bm{\mathcal{A}}`; composed propagators reproduce a single
  ``expm``; the uniform-field shortcut agrees with the stepped path; the
  photon collapse switches exactly at :math:`\kappa d = 12`;
* The command line — every subcommand parses, runs, and writes a file whose
  header describes its own columns.

**Invariants and regressions.**  Pinned values are the weakest kind of test:
a wrong curve that is merely self-consistent will pass them.  Where possible
the suite asserts a contract instead — for example that branch 1 is the
*globally lowest* root at every :math:`pd`, which is what the warm start
once violated.  A small number of reference values guard against silent
numerical drift.
