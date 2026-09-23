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

Paschen's law [Raizer1991]_ states that the DC breakdown voltage of a uniform-field gap is
a universal function of the product :math:`pd` of gas pressure and electrode
separation.  It follows from Townsend's theory: avalanche growth is
parameterised by the first Townsend coefficient :math:`\alpha(E/N)` together
with a cathode secondary-electron-emission (SEE) efficiency :math:`\gamma`,
giving the classical criterion :math:`\alpha d = \ln(1 + \gamma^{-1})`.

Paschen's law is known to be inaccurate beyond
:math:`pd \gtrsim 1\text{–}10` bar·mm.  In electronegative gases such as air
part of that inaccuracy comes from a mechanism that the classical theory
omits: **negative-ion transit and detachment**.  Electron attachment does not
necessarily remove an electron permanently.  If a negative ion drifts across
the gap and only *later* detaches (collisionally or associatively), it
returns its electron to the swarm mid-gap, where it keeps contributing to
avalanche growth.  Whether this matters depends on how the negative-ion
transit time :math:`d/v_-` compares with the detachment time — a question
that a purely local ionization coefficient :math:`\alpha(E/N) - \eta(E/N)`
cannot answer.

``Incept1D`` therefore treats the gap as a boundary-value problem for the
coupled fluxes of electrons, positive ions, several negative-ion species,
and ionizing photons, and reduces "does a self-sustained discharge start in
this gap?" to the vanishing of a determinant.

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
computed by the same routine, :func:`incept1d.solver.inception_det`.  The
derivation is spread over the following pages:

* :ref:`Chap:Transport` — the augmented ODE and its formal solution;
* :ref:`Chap:Photoionization` — the two-stream photon equations;
* :ref:`Chap:SecondaryEmission` — the cathode boundary condition;
* :ref:`Chap:InceptionCriterion` — assembling :math:`\bm{Q}` and recovering
  Paschen's law;
* :ref:`Chap:AirScheme` — the dry-air reaction scheme shipped with the code.
