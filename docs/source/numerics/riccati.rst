.. _Chap:Numerics:Riccati:

The reflection (Riccati) criterion
==================================

.. contents:: On this page
   :local:
   :depth: 1

The inception condition :eq:`eq_det_criterion` can be evaluated without
ever forming the propagator :math:`\bm{M}(d)`.  This is the default
(``--criterion riccati``, :func:`incept1d.solver.riccati_criterion`); the
determinant itself (``--criterion detq``, :ref:`Chap:Numerics:Determinant`)
remains available and solves the same discrete problem on the same grid.

Why not form the propagator
---------------------------

:eq:`eq_Qd` imposes the anode conditions on :math:`\bm{M}(d)\vec\theta_0`,
i.e. on the state carried from the cathode across the whole gap.  Half of
the state, though, travels the other way: positive ions drift towards the
cathode, and the backward photon streams :math:`\Psi^-_j` have their
boundary condition at the anode.  Integrated from the cathode, these
components grow like :math:`e^{+\kappa_j x}` and, at :math:`\lambda > 0`,
like :math:`e^{+\lambda x/|v_+|}`.  Once the modes the anode rows need
differ in growth by more than about 35 e-folds, every row of
:math:`\bm{Q}_d` is parallel to the dominant mode to machine precision,
and :math:`\det\bm{Q}` can no longer be evaluated from :math:`\bm{M}(d)`.
An optically thick photon group does that on its own: :math:`\kappa_j d`
of several hundred is common.  The difficulty is the direction of
integration, not the physics, and it disappears when every component is
integrated in the direction it travels.

The reflection operator
-----------------------

Split the augmented state into the components whose conditions are
imposed at the cathode, :math:`\vec\theta_f` (electrons, negative ions,
forward photons), and those whose conditions are imposed at the anode,
:math:`\vec\theta_b` (positive ions, backward photons: the components
:math:`\bm{\Pi}_+` and :math:`\bm{\Pi}_{\Psi^-}` select).  With
:math:`\bm{\mathcal{A}}` partitioned accordingly, the anode condition
:math:`\vec\theta_b(d) = \vec 0` is carried back to the cathode as

.. math::
   :label: eq_reflection

   \vec\theta_b(x) = \bm{P}(x)\,\vec\theta_f(x), \qquad \bm{P}(d) = \bm{0},

where :math:`\bm{P}` obeys the matrix Riccati equation

.. math::
   :label: eq_riccati

   \bm{P}' = \bm{\mathcal{A}}_{bf} + \bm{\mathcal{A}}_{bb}\bm{P}
             - \bm{P}\bm{\mathcal{A}}_{ff} - \bm{P}\bm{\mathcal{A}}_{fb}\bm{P},

integrated from :math:`x = d` down to :math:`x = 0`.  In that direction
every component travels the way it physically does, so the modes that
make :math:`\bm{M}` useless decay instead of growing.  :math:`\bm{P}` is
the reflection operator of the part of the gap downstream of :math:`x`:
the fluxes that come back towards the cathode per unit flux going
towards the anode.

At the cathode, :math:`\bm{Q}_0\vec\theta_0 = \vec 0` with
:math:`\vec\theta_b = \bm{P}(0)\vec\theta_f` leaves a small square
system, and the inception condition becomes

.. math::
   :label: eq_riccati_criterion

   g = \det\bigl(\bm{Q}_{0,f} + \bm{Q}_{0,b}\,\bm{P}(0)\bigr) = 0 .

:math:`g` vanishes exactly where :math:`\det\bm{Q}` does: both state that
the cathode and anode conditions admit a common non-trivial solution.
With a single electron-emission row, :math:`g` is one minus the loop gain
— one minus the secondary electrons produced per electron leaving the
cathode, through every channel, with :math:`\bm{P}(0)\vec e_e` the ion
and backward-photon fluxes that return.  It is therefore of order one,
positive below inception and negative above for every mechanism, whereas
the sign of :math:`\det\bm{Q}` depends on the ordering of its rows.

Solving one step exactly
------------------------

On each integration step :math:`\bm{\mathcal{A}}` is held constant (the
same midpoint or Magnus exponent :math:`\bm\Omega` that builds
:math:`\bm{M}`, :ref:`Chap:Numerics:Propagator`), and
:eq:`eq_riccati` is then solved exactly.  With
:math:`\bm{E} = \exp(-\bm\Omega)` relating the two ends of the step, the
step is described by its reflection and transmission matrices,

.. math::

   \bm{T}_f = \bm{E}_{ff}^{-1}, \quad
   \bm{R}_r = -\bm{E}_{ff}^{-1}\bm{E}_{fb}, \quad
   \bm{R}_l = \bm{E}_{bf}\bm{E}_{ff}^{-1}, \quad
   \bm{T}_b = \bm{E}_{bb} - \bm{E}_{bf}\bm{E}_{ff}^{-1}\bm{E}_{fb},

and the reflection operator is carried across it by

.. math::
   :label: eq_reflection_step

   \bm{P} \leftarrow \bm{R}_l + \bm{T}_b\,\bm{P}\,
                     (\bm{I} - \bm{R}_r\bm{P})^{-1}\,\bm{T}_f .

These four matrices are formed directly only for a thin substep whose
propagator has a condition number below :math:`10^{12}`.  A whole step is
built from such substeps by the adding–doubling method
(:func:`incept1d.solver._slab_scattering`): two identical slabs combine by
the Redheffer star product (:func:`incept1d.solver._star`), so a step of
:math:`m` substeps takes :math:`\log_2 m` combinations.  That matters
where the drift of a slow species is small and its reaction rates, divided
by the drift speed, span millions of e-folds per step.  For a system whose
couplings are all non-negative sources the reflection and transmission
matrices are non-negative too, and composing them involves no
cancellation — the property that makes this stable where re-orthonormalising
the propagated subspace is not (a transiently tiny component, such as
electrons attaching in a low-field region before an avalanche downstream,
would be buried under the rounding error of the large ones).

A step whose transmissions overflow has an astronomically large gain, and
is read as above inception.

Poles
-----

:math:`\bm{P}` stays finite below inception.  It can only blow up where a
part of the gap, :math:`[x, d]`, sustains itself without the cathode —
through photoionization alone, for example — and then the whole gap is
supercritical.  That shows as the loop gain :math:`\bm{R}_r\bm{P}` of
:eq:`eq_reflection_step` (or :math:`\bm{R}_r^A\bm{R}_l^B` inside a star
product) reaching spectral radius one, and the criterion returns
:math:`-1`.  The test is on the spectral radius, not on a sign change of
:math:`\det(\bm{I} - \bm{R}_r\bm{P})`: for non-negative loop gains
:math:`\bm{I} - \bm{K}` has a non-negative inverse exactly when
:math:`\rho(\bm{K}) < 1`, and :math:`\rho` stays above one once it has
crossed, whereas the determinant changes sign back when a second mode
crosses too — which far above inception happens within one step.

Inception can therefore show in two ways: as a zero of :math:`g`, or,
when the cathode contributes little, as the pole at which :math:`g` jumps
from about :math:`+1` to :math:`-1`.  Both are crossings from :math:`+` to
:math:`-` in increasing :math:`E/N`, and the root finders accept exactly
those (:ref:`Chap:Numerics:RootFinding`).  In :math:`\lambda`, :math:`g`
rises monotonically — every loop gain falls as :math:`\lambda` grows — so
its single sign change is the dominant mode.

Verification
------------

``--verify`` evaluates :math:`\det\bm{Q}` on the same field and grid just
below and just above every root found and reports whether it changes
sign.  The test suite checks the criterion against the exact
:math:`\det\bm{Q}` of a uniform gap computed in high-precision arithmetic,
against :math:`\det\bm{Q}` on matched grids in non-uniform gaps, and
against an independent discretisation of the time-dependent equations
whose fastest mode must equal the growth rate :math:`\lambda^*`.

Cost
----

On the same grid the criterion costs about as much as forming
:math:`\bm{M}` and :math:`\det\bm{Q}` directly: one small matrix
exponential per step, plus a few operations on matrices whose size is the
number of components of each kind.  It does not grow with the optical
depth of the photon groups, which is what makes the determinant route
expensive once every group is explicit.
