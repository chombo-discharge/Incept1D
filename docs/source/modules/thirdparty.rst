.. _Chap:ThirdParty:

Third-party integrations
========================

.. contents:: On this page
   :local:
   :depth: 1

A mechanism file is a complete description of a gas chemistry, and the 1-D
criterion is rarely the last word: once a gap is known to break down, the
next question is usually what the discharge does in three dimensions.  The
commands on this page export a mechanism in the form an external code
expects, so that the *same* chemistry drives both calculations.

Only ``chombo-discharge`` is supported today.  Each further integration
gets its own subcommand and its own section here rather than changing the
mechanism interface, which stays independent of any particular consumer.

.. _Chap:CreateChomboDischargeData:

``chombo-discharge`` — ``incept1d chombo``
------------------------------------------

`chombo-discharge <https://github.com/chombo-discharge/chombo-discharge>`_
[Marskar2019]_ is a 3-D plasma fluid solver.  ``incept1d chombo`` exports
the transport and rate-coefficient tables it needs as a function of
:math:`E/N`.  It does not solve the inception problem itself.

Inputs
~~~~~~

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
~~~~~~~

.. code-block:: bash

   incept1d chombo mechanisms/air/pancheshnyi/air_pancheshnyi.py \
       --modifier mechanisms/air/pancheshnyi/databases.json --config-label Phelps \
       --min-EN 1 --max-EN 2000 --num-EN 500 --write-to-file air_phelps.dat

.. note::

   This command has no ``--no-plot`` flag; on headless machines set
   ``MPLBACKEND=Agg``.

Outputs
~~~~~~~

A multi-panel sanity-check figure — ionization and attachment coefficients,
mean energy, mobilities and diffusion coefficients, and every rate
coefficient against :math:`E/N` — so that a new or modified mechanism can be
inspected before it is handed to the 3-D solver.

``--write-to-file`` writes the table itself: one row per :math:`E/N` point,
with a self-describing header naming the species order, the reaction order
and the units of every column:

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

Everything is reduced by the number density (``alpha/N``, ``mu*N``,
``D*N``), so one table serves every pressure the 3-D solver runs at.

API reference
~~~~~~~~~~~~~

.. automodule:: incept1d.chombo
   :members:
   :undoc-members:
