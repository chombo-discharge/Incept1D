.. _Chap:Eigenvalues:

Local growth modes — ``incept1d eigenvalues``
=============================================

.. contents:: On this page
   :local:
   :depth: 1

``incept1d eigenvalues`` computes the eigenvalues of the *local* transport matrix
:math:`\bm{A} = \bm{R}\bm{V}^{-1}` as a function of :math:`E/N` at a single
spatial point — no gap integration, no boundary conditions.  It is the
purely local counterpart of the full inception criterion: a positive real
eigenvalue means that, locally, net ionization exceeds attachment/loss and
the charge-carrier flux would grow exponentially in an infinite uniform
medium.  As shown in :ref:`Chap:InceptionCriterion`, the largest eigenvalue
reduces to :math:`\lambda_+` of the closed-form three-species model and is
the *apparent effective ionization coefficient*.

Eigenvalues have no inherent ordering, so they are matched between
adjacent :math:`E/N` samples before plotting: a track follows one mode
across the sweep instead of jumping to another wherever two of them cross.

Inputs
------

.. code-block:: console

   incept1d eigenvalues MECHANISM.py [CONFIG.json ...]
       [--p P] [--T T] [--EN-lo LO] [--EN-hi HI] [--EN-num N]
       [--pressure-scan] [--p-min MIN] [--p-max MAX] [--p-num N] [--eig-index IDX]
       [--write-to-file FILE]

Two modes are available.  The **default mode** plots every eigenvalue track
(:math:`\mathrm{Re}\,\lambda_j/N` vs. :math:`E/N`) at a fixed pressure, one
subplot per configuration.  The **pressure-scan mode** (``--pressure-scan``)
plots one selected track (``--eig-index``, default 0 = leading mode) vs.
:math:`E/N` for several log-spaced pressures between ``--p-min`` and
``--p-max``.  Because the three-body attachment and the detachment reactions
scale differently with :math:`N`, the field at which the leading eigenvalue
crosses zero is pressure dependent — the pressure scan shows this directly.

Examples
--------

.. code-block:: bash

   # All eigenvalue tracks for the baseline and the no-detachment variant
   incept1d eigenvalues mechanisms/air/air_pancheshnyi.py mechanisms/air/nodetachment.json --EN-lo 20 --EN-hi 300

   # Leading eigenvalue vs E/N for five pressures between 1 mbar and 10 bar
   incept1d eigenvalues mechanisms/air/air_pancheshnyi.py --pressure-scan --p-min 1e-3 --p-max 10 --p-num 5

.. note::

   ``incept1d eigenvalues`` has no ``--no-plot`` flag; on headless machines set
   ``MPLBACKEND=Agg`` and use ``--write-to-file``.

Outputs
-------

A figure with one subplot per configuration: the reduced eigenvalues
:math:`\mathrm{Re}(\lambda_j/N)` in m\ :sup:`2` against :math:`E/N`, one
line per track, or one line per pressure in ``--pressure-scan`` mode.
Dividing by the number density is what makes the curves collapse across
pressures for a chemistry whose rates scale with :math:`N`.

``--write-to-file`` writes one tab-separated file **per configuration**,
named after the label, with a metadata header and columns ``EN_Td``
followed by ``Re_lam<j>_N_m2`` for each track.

API reference
-------------

.. automodule:: incept1d.eigenvalues
   :members:
   :undoc-members:
