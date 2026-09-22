.. _Chap:CreateChomboDischargeData:

CreateChomboDischargeData.py — transport tables
===============================================

.. contents:: On this page
   :local:
   :depth: 1

``CreateChomboDischargeData.py`` is the bridge between a mechanism file and
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

.. literalinclude:: ../../../CreateChomboDischargeData.py
   :language: python
   :pyobject: generate

Command-line interface
----------------------

.. code-block:: console

   python3 CreateChomboDischargeData.py MECHANISM.py
       [--modifier CONFIG.json] [--config-label LABEL]
       [--min-EN TD] [--max-EN TD] [--num-EN N]
       [--pressure BAR] [--temperature K] [--write-to-file FILE]

Unlike the solver scripts, this script takes a single configuration through
``--modifier`` (a JSON file in the same format as :ref:`Chap:Configuration`)
and ``--config-label`` (which entry to use).  The :math:`E/N` grid defaults
to 1–10\ :sup:`6` Td with 1000 log-spaced points.

Example
-------

.. code-block:: bash

   python3 CreateChomboDischargeData.py Air/Air_Pancheshnyi.py \
       --modifier Air/Databases.json --config-label Phelps \
       --min-EN 1 --max-EN 2000 --num-EN 500 --write-to-file air_phelps.dat

.. note::

   This script has no ``--no-plot`` flag; on headless machines set
   ``MPLBACKEND=Agg``.

API reference
-------------

.. automodule:: CreateChomboDischargeData
   :members:
   :undoc-members:
