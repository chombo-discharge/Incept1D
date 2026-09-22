.. _Chap:CreateChomboDischargeData:

Transport tables — ``incept1d chombo``
======================================

.. contents:: On this page
   :local:
   :depth: 1

``incept1d chombo`` is the bridge between a mechanism file and
the external 3-D plasma solver
`chombo-discharge <https://github.com/chombo-discharge/chombo-discharge>`_
[Marskar2019]_.  It does not solve the inception problem — it exports the
raw transport and rate-coefficient tables a fluid solver needs as a
function of :math:`E/N`, so that the *same* chemistry used for the 1-D
criterion can be used in a 3-D simulation of the onset.

Output columns
--------------

.. code-block:: text

   1      E/N [Townsend]
   2      alpha/N [m^2]
   3      eta/N [m^2]
   4      Electron mean energy [eV]
   5-6    mu*N, D*N for e
   7-8    mu*N, D*N for the first ion species
   ...    (two columns per species, in SPECIES order)
   next   Raw rate coefficient for each reaction, in REACTIONS order
          (m^3/s for two-body, m^6/s for three-body; flagged in the header)

:func:`generate` evaluates the mechanism over an :math:`E/N` grid;
:func:`build_header` writes a self-describing header (species order,
reaction order, units); :func:`plot_results` renders a multi-panel
sanity-check figure (ionization/attachment coefficients, mean energy,
mobilities and diffusion coefficients, rate coefficients) so a new or
modified mechanism can be inspected before it is handed to the 3-D solver.

.. literalinclude:: ../../../src/incept1d/chombo.py
   :language: python
   :pyobject: generate

Command-line interface
----------------------

.. code-block:: console

   incept1d chombo MECHANISM.py
       [--modifier CONFIG.json] [--config-label LABEL]
       [--min-EN TD] [--max-EN TD] [--num-EN N]
       [--pressure BAR] [--temperature K] [--write-to-file FILE]

Unlike the solver commands, ``incept1d chombo`` takes a single configuration through
``--modifier`` (a JSON file in the same format as :ref:`Chap:Configuration`)
and ``--config-label`` (which entry to use).  The :math:`E/N` grid defaults
to 1–10\ :sup:`6` Td with 1000 log-spaced points.

Example
-------

.. code-block:: bash

   incept1d chombo mechanisms/air/air_pancheshnyi.py \
       --modifier mechanisms/air/databases.json --config-label Phelps \
       --min-EN 1 --max-EN 2000 --num-EN 500 --write-to-file air_phelps.dat

.. note::

   This command has no ``--no-plot`` flag; on headless machines set
   ``MPLBACKEND=Agg``.

API reference
-------------

.. automodule:: incept1d.chombo
   :members:
   :undoc-members:
