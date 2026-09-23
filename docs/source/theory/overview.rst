.. _Chap:TheoryOverview:

The one-dimensional model
=========================

.. note::

   The theory pages derive the model that the code implements, from the
   governing equations to the determinant inception criterion, and map
   every equation onto the function that implements it.  Literature
   sources for the physical data and approximations are listed in
   :ref:`Chap:References`.

.. contents:: On this page
   :local:
   :depth: 1

Motivation
----------

**Paschen's law** [Raizer1991]_ states that the DC breakdown voltage of a
uniform-field gap is a universal function of the product :math:`pd` of gas
pressure and electrode separation.  It follows from Townsend's theory:
avalanche growth is parameterised by the first Townsend coefficient
:math:`\alpha(E/N)` together with a cathode secondary-electron-emission (SEE)
efficiency :math:`\gamma`, giving the classical criterion
:math:`\alpha d = \ln(1 + \gamma^{-1})`.  The law is exact within its own
assumptions, and :ref:`Chap:Examples:Paschen` reproduces it to solver
tolerance.  Those assumptions are restrictive: a uniform field, a single
ionizing process, and a discharge sustained entirely by feedback at the
cathode.

**The streamer criterion** [Raizer1991]_ addresses the opposite regime.  At
larger :math:`pd` a single avalanche can accumulate enough space charge to
distort the applied field and propagate on its own, without waiting for the
cathode.  Breakdown is then declared when the avalanche reaches a critical
size, :math:`\int\max(\alpha - \eta, 0)\,dx \approx 18`\ –\ :math:`20`.
This is a local criterion evaluated along the field line: it asks only how
much net ionization an avalanche accumulates, not where the charge ends up or
what becomes of it, and the threshold constant is empirical.

``Incept1D`` avoids choosing between them.  Rather than asking whether an
avalanche reaches a given size, it treats the gap as a boundary-value problem
for the coupled fluxes of electrons, positive ions, several negative-ion
species and ionizing photons, and asks whether a self-sustaining mode exists
at all.  Negative-ion transit, detachment and ion conversion, two-stream
photoionization, and both ion- and photon-induced secondary emission are
carried as part of that system rather than folded into an effective
coefficient, so cathode feedback and volume feedback are weighed against each
other instead of being assumed.  The question "does a self-sustained discharge
start in this gap?" reduces to the vanishing of a determinant.

Species and governing equation
------------------------------

The full three-dimensional evolution of electrons and ions is approximated
by a one-dimensional drift-reaction equation along the field line:

.. math::
   :label: eq_drift_reaction

   \partial_t\vec{n} + \partial_x\left(\bm{V}\vec{n}\right)
     = \bm{R}\vec{n} + \sum_j\vec{\beta}_j\,\xi_j\,\kappa_j\,\Psi_j^0 .

Here

* :math:`\vec{n}(x,t)` is the column vector of the :math:`N_s` tracked
  species number densities;
* :math:`\bm{V} = \diag(v_1, v_2, \ldots)` is the diagonal matrix of signed
  drift velocities;
* :math:`\bm{R}` is the reaction-rate matrix, which crucially does **not**
  depend on :math:`\vec{n}`;
* The last term is photoionization by :math:`N_\gamma` photon groups:
  :math:`\kappa_j\Psi_j^0` is the number of group-:math:`j` photons absorbed
  per unit volume and time, :math:`\xi_j` the photoionization efficiency,
  and :math:`\vec{\beta}_j` selects which species the absorption produces.

Two approximations are built in at this level:

**Local field approximation (LFA).**  The rates in :math:`\bm{R}` and the
velocities in :math:`\bm{V}` depend only on the local reduced field
:math:`E/N`, the gas pressure, and the temperature.  Electron transport and
rate coefficients are therefore taken from a Boltzmann solver (BOLSIG+)
tabulated against :math:`E/N`.

**Linearity.**  Because :math:`\bm{R}` is independent of :math:`\vec{n}`,
equation :eq:`eq_drift_reaction` is linear and homogeneous.  Space-charge
effects, electron-ion recombination, and any other density-dependent
processes are excluded.  This is appropriate for the *inception* stage, where
densities are small, and is what makes the eigenvalue reduction in
:ref:`Chap:Transport` possible.

What the model neglects
-----------------------

* **Charged-species diffusion** is ignored [Ferreira2019]_.  Diffusion is mostly relevant
  where strong gradients exist (sheaths, wall losses, back-diffusion to the
  cathode), not in the bulk of a growing discharge.  Including it is
  possible but requires a more sophisticated numerical treatment.
* **Transverse structure.**  The model lives on a single field line.  For
  non-uniform fields (sphere-plane, sphere-sphere, tabulated field lines)
  the on-axis or along-line field is used and flux-tube divergence is
  neglected; see :ref:`Chap:FieldDistributions`.
* **Off-axis photons.**  Radiative transfer is solved in a two-stream
  approximation along :math:`\pm x` only, with a free cone-angle parameter
  that compensates for the fraction of photons emitted in other directions;
  see :ref:`Chap:Photoionization`.

What the model predicts
-----------------------

Seeking solutions of the form :math:`\vec{n}(x,t) = \vec{n}(x)e^{\lambda t}`
turns :eq:`eq_drift_reaction` into a linear ODE in :math:`x`, whose formal
solution is a (path-ordered) matrix exponential.  Imposing the cathode and
anode boundary conditions — zero incoming ion fluxes, secondary emission at
the cathode, no incoming photons — gives a homogeneous linear system
:math:`\bm{Q}(\lambda)\vec{\theta}_0 = \vec{0}` for the cathode fluxes
:math:`\vec{\theta}_0`.  A non-trivial solution, i.e. a discharge mode,
exists if and only if

.. math::

   \det\bm{Q}(\lambda) = 0 .

The **inception threshold** is the applied field at which this holds for
:math:`\lambda = 0`; the **temporal growth rate** above threshold is the
:math:`\lambda > 0` at which it holds for a fixed applied field.  Both are
computed by the same routine, :func:`incept1d.solver.inception_det`.

Where the derivation continues
------------------------------

The rest of this chapter builds :math:`\bm{Q}` one piece at a time:

* :ref:`Chap:Transport` — The charged-species equation, in flux form.
* :ref:`Chap:Photoionization` — The two-stream photon equations.
* :ref:`Chap:SecondaryEmission` — What happens at the two electrodes.
* :ref:`Chap:AugmentedSystem` — The three combined into one ODE, and its
  formal solution.
* :ref:`Chap:InceptionCriterion` — Assembling :math:`\bm{Q}`, and recovering
  Paschen's law from it as a limiting case.

The reaction scheme itself is not part of the derivation: it is an input.
The dry-air scheme shipped with the code is documented as a worked example in
:ref:`Chap:AirScheme`.

How the assembled criterion is then evaluated and its roots located is a
numerical question, taken up in :ref:`Chap:Numerics`.
