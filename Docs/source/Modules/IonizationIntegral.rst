.. _Chap:IonizationIntegral:

Ionization integrals — ``incept1d ionization``
==============================================

.. contents:: On this page
   :local:
   :depth: 1

``incept1d ionization`` evaluates the two local surrogates of the
inception criterion defined in :ref:`Chap:InceptionCriterion` along the
(possibly non-uniform) field profile and plots them against applied
voltage:

* the effective ionization integral :math:`I_\alpha = \int_0^d
  \max(\alpha - \eta, 0)\,dx` (:eq:`eq_ionization_integral`), i.e. the
  classical streamer-criterion integrand, which ignores detachment, ion
  transit and photoionization entirely;
* the apparent effective ionization integral :math:`I_\lambda = \int_0^d
  \max(\mathrm{Re}\,\lambda_{\max}(\bm{R}\bm{V}^{-1}), 0)\,dx`
  (:eq:`eq_apparent_ionization_integral`), which uses the leading local
  eigenvalue in place of :math:`\alpha - \eta`.

Comparing the two curves at the inception voltage found by ``incept1d pdiv``
quantifies how much the local approximation misses once negative-ion
transit and detachment matter, and how sensitive that is to the choice of
electron cross sections.

.. literalinclude:: ../../../src/incept1d/ionization.py
   :language: python
   :pyobject: aed_integral

.. literalinclude:: ../../../src/incept1d/ionization.py
   :language: python
   :pyobject: eig_integral

Command-line interface
----------------------

.. code-block:: console

   incept1d ionization MECHANISM.py [CONFIG.json ...]
       --d D [D ...] (--p P [P ...] | --data-file FILE [--pressure-column COL] [--voltage-column COL])
       [--voltage-lo V] [--voltage-hi V] [--voltage-num N] [--single-voltage V]
       [--T T] [--field SPEC] [--N N] [--no-plot] [--write-to-file FILE]

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Meaning
   * - ``--d D [D ...]``
     - Gap distance(s) in mm (required).
   * - ``--p P [P ...]``
     - Pressure(s) in bar (required unless ``--data-file``).
   * - ``--voltage-lo/-hi/-num``
     - Linear voltage sweep in kV (defaults 0–1 kV, 100 points).
   * - ``--single-voltage V``
     - Evaluate at one voltage, print the result, no plot.
   * - ``--data-file FILE``
     - Read (pressure, voltage) pairs from an ASCII table — e.g. measured
       breakdown voltages — and evaluate the integrals at exactly those
       points instead of sweeping.  Columns are selected by name or 0-based
       index.
   * - ``--N N``
     - Number of midpoint-rule quadrature steps across the gap for
       non-uniform fields (default 200).
   * - ``--field SPEC``
     - Geometry, as for ``incept1d pdiv``.

Examples
--------

.. code-block:: bash

   # Both integrals vs voltage, four cross-section sets, 10 mm gap at 10 bar
   incept1d ionization mechanisms/Air/Air_Pancheshnyi.py mechanisms/Air/Databases.json \
       --p 10 --d 10 --voltage-lo 100 --voltage-hi 300

   # Evaluate the integrals at measured breakdown voltages in a sphere-plane gap
   incept1d ionization mechanisms/Air/Air_Pancheshnyi.py --d 20 \
       --data-file measurements.dat --pressure-column 0 --voltage-column 1 \
       --field sphere-plane 25

API reference
-------------

.. automodule:: incept1d.ionization
   :members:
   :undoc-members:
