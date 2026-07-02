Inception
=========

``Inception.py`` is the core solver: it implements the boundary-value
inception criterion :math:`\det\bm{Q}(\lambda=0) = 0` derived in
:doc:`../model` and provides the CLI used to compute generalized Paschen
curves.

Loading a mechanism
--------------------

:func:`load_mechanism` reads a mechanism ``.py`` file (e.g.
``Air/Air_Pancheshnyi.py``) with ``importlib`` — mechanism files are *not*
imported as ordinary packages, they are executed as standalone modules, with
optional override variables injected into their namespace beforehand (via a
companion ``Config.py`` living next to the mechanism file; see
``Air/Config.py`` for the protocol it implements). The loaded module is
validated against a fixed required interface and wrapped in a
:class:`Mechanism` object, which bakes in every configuration parameter
(photoionization/photoemission scale factors, reaction-rate multipliers,
per-polarity SEE overrides) so that every other function in the repository
can treat ``mod.get_R(EN, p, T)`` etc. as the final, fully-configured
answer.

.. literalinclude:: ../../Inception.py
   :language: python
   :pyobject: load_mechanism

.. autoclass:: Inception.Mechanism
   :members:

Assembling the augmented ODE
------------------------------

:func:`_build_A_aug` builds the block matrix :math:`\bm{\mathcal{A}}` from
``eq:augmented_ode`` (see :doc:`../model`) at a given reduced field, gas
pressure/temperature, and temporal growth rate :math:`\lambda`. Photon
groups whose optical depth over the gap is very large are folded into a
*local* absorption term added directly to :math:`\bm{A}` (a numerical
convenience — those groups are re-absorbed within a fraction of a step, so
propagating them explicitly would only add stiffness, not physics).

.. literalinclude:: ../../Inception.py
   :language: python
   :pyobject: _build_A_aug

Boundary conditions and the determinant
------------------------------------------

:func:`_assemble_det_Q` builds :math:`\bm{Q} = [\bm{Q}_0; \bm{Q}_d]` from the
propagator matrix :math:`\bm{M}(d)` and the mechanism's row-selection /
SEE-efficiency functions, and returns :math:`\det\bm{Q}` (as a signed,
overflow-guarded log-determinant; ``NaN`` signals an ill-conditioned or
overflowing system rather than a numerically meaningless value).

.. literalinclude:: ../../Inception.py
   :language: python
   :pyobject: _assemble_det_Q

:func:`inception_det` ties the two together: it integrates
:math:`\bm{\mathcal{A}}(x)` across the gap (adaptively refining the grid
when using :func:`midpoint_propagator`) to build :math:`\bm{M}(d)`, then
calls ``_assemble_det_Q``. This is the function whose root in :math:`E/N`,
at fixed :math:`p \cdot d`, *is* the inception voltage.

Finding roots and tracking Paschen branches
-----------------------------------------------

:func:`find_all_breakdown_EN` brackets and root-finds all :math:`E/N`
solutions of ``inception_det(...) == 0`` at a single :math:`p\cdot d`, and
:func:`compute_paschen_curve` sweeps this over an array of :math:`p\cdot d`
values, tracking each solution as a continuous *branch* (matching new roots
to existing branches by nearest-neighbour distance in :math:`\log(E/N)`, so
that saddle-node bifurcations that create new branches mid-sweep don't
corrupt existing ones).

.. literalinclude:: ../../Inception.py
   :language: python
   :pyobject: compute_paschen_curve

Command-line usage
--------------------

.. code-block:: console

   python Inception.py <mechanism.py> [CONFIG.json ...] \
       [--p P [P ...]] [--d D [D ...]] [--T T] \
       [--pd-min PD] [--pd-max PD] [--pd-num N] \
       [--field SPEC] [--xi-photo tag=val] [--xi-emit tag=val]

See the module docstring (below) for the full argument reference.

API reference
--------------

.. automodule:: Inception
   :members:
   :undoc-members:
   :exclude-members: Mechanism
