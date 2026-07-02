CreateChomboDischargeData
==========================

``CreateChomboDischargeData.py`` is the bridge between a mechanism file and
the external 3-D plasma solver ``chombo-discharge`` (the "three-dimensional
plasma simulations" used in the manuscript to verify the spatio-temporal
onset predicted by the 1-D model). It does not solve the inception problem
itself — it exports the raw transport and rate-coefficient tables a 3-D
solver needs as a function of :math:`E/N`:

.. code-block:: text

   1      E/N [Townsend]
   2      alpha/N [m^2]
   3      eta/N [m^2]
   4      Electron mean energy [eV]
   5-6    mu*N, D*N for e
   ...    mu*N, D*N for each ion species, in SPECIES order
   next   Raw rate coefficient for each reaction in REACTIONS order

:func:`generate` evaluates the mechanism's rates and transport coefficients
over an :math:`E/N` grid; :func:`build_header` writes a self-describing
column header (species order, reaction order, three-body vs. two-body
units); :func:`plot_results` renders a four-panel sanity-check figure
(ionization/attachment coefficients, mean energy, mobilities/diffusion
coefficients, rate coefficients) so a new or modified mechanism can be
visually checked before being handed to the 3-D solver.

.. literalinclude:: ../../CreateChomboDischargeData.py
   :language: python
   :pyobject: generate

Command-line usage
--------------------

.. code-block:: console

   python CreateChomboDischargeData.py MECHANISM_FILE \
       [--modifier MODIFIER_JSON] [--config-label LABEL] \
       [--min-EN 10] [--max-EN 10000] [--num-EN 500] \
       [--pressure 1.0] [--temperature 300.0] \
       [--output OUTPUT_FILE]

API reference
--------------

.. automodule:: CreateChomboDischargeData
   :members:
   :undoc-members:
