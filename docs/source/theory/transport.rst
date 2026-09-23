.. _Chap:Transport:

The drift-reaction system
=========================

.. contents:: On this page
   :local:
   :depth: 1

This page reduces the governing equation :eq:`eq_drift_reaction` to a
first-order ODE in :math:`x` for the charged-species fluxes.  The photon
fluxes that appear as a source term are derived in
:ref:`Chap:Photoionization`, the electrode conditions in
:ref:`Chap:SecondaryEmission`, and the three are assembled into one system in
:ref:`Chap:AugmentedSystem`.

Coordinate convention
---------------------

The spatial coordinate :math:`x` runs from the **cathode at** :math:`x = 0`
**to the anode at** :math:`x = d`.  The electric field points in the
:math:`-x` direction, but the code works exclusively with the magnitude
:math:`E/N > 0`.

The direction each species travels is therefore carried by the sign of its
entry in the diagonal drift-velocity matrix :math:`\bm{V}`, which the
mechanism supplies directly:

.. math::

   V_{ii} > 0 \quad\text{(electrons and negative ions, towards the anode)},
   \qquad
   V_{ii} < 0 \quad\text{(positive ions, towards the cathode)}.

Nothing in the solver assumes how those velocities were obtained.  A
mechanism may return a mobility times the field, a tabulated drift velocity
from a Boltzmann solver, or a fitted expression — only the signed velocity
itself is part of the interface (:ref:`Chap:NewMechanisms`).  With this
convention a positive eigenvalue of :math:`\bm{R}\bm{V}^{-1}` corresponds to
a spatially growing solution; see :ref:`Chap:Eigenvalues`.

For a non-uniform gap the field varies along the line as
:math:`E(x) = E_\mathrm{ref}\,f(x/d)`, where :math:`f` is the normalised
geometry profile defined in :ref:`Chap:FieldDistributions`.  Everything below
holds pointwise in :math:`x`; the only consequence is that
:math:`\bm{R}` and :math:`\bm{V}`, which depend on the local :math:`E/N`,
become functions of position.

Reduction to an ODE
-------------------

Starting from :eq:`eq_drift_reaction`, follow the eigenvalue approach of
[Ferreira2019]_ and seek solutions with exponential time dependence,
:math:`\vec{n}(x,t) \to \vec{n}(x)e^{\lambda t}` and
:math:`\vec{\Psi}(x,t) \to \vec{\Psi}(x)e^{\lambda t}`.  This is allowed
without loss of generality because the system is linear and homogeneous.

Writing the result in terms of the **charged-species flux**
:math:`\vec{y} = \bm{V}\vec{n}` rather than the density gives

.. math::
   :label: eq_flux_ode

   \partial_x\vec{y} = \left(\bm{R} - \lambda\bm{I}\right)\bm{V}^{-1}\vec{y}
     + \sum_j\vec{\beta}_j\xi_j\kappa_j\left(\Psi_j^+ + \Psi_j^-\right),

where :math:`\Psi_j^\pm` are the forward and backward photon fluxes of group
:math:`j`.  They are not yet determined — :eq:`eq_flux_ode` is only half of
the system, and :ref:`Chap:Photoionization` supplies the other half.  A
mechanism with no photoionization drops the sum entirely.

Why fluxes and not densities
............................

Working with :math:`\vec{y}` rather than :math:`\vec{n}` is what keeps the
problem linear *and* the boundary conditions simple.  Every condition at an
electrode is a statement about what crosses it: no positive ions arrive at
the anode, no negative ions arrive at the cathode, the electron flux leaving
the cathode is whatever secondary emission produces.  Each is a linear
constraint on :math:`\vec{y}` evaluated at a single point, which is exactly
what a boundary-value problem needs.  In terms of densities the same
statements would carry the drift velocities around with them.

The transport matrix and its eigenvalues
----------------------------------------

Set :math:`\lambda = 0` and drop the photon term.  What remains is

.. math::

   \partial_x\vec{y} = \bm{R}\bm{V}^{-1}\vec{y},

so the eigenvalues :math:`k_i` of :math:`\bm{R}\bm{V}^{-1}` are the *local*
spatial growth rates of the flux modes: the solution is a sum of eigenmodes
:math:`\sum_i c_i\vec{l}_ie^{k_ix}`, and if any :math:`\mathrm{Re}\,k_i > 0`
there is net spatial growth.

The largest such eigenvalue, :math:`\lambda_{\max}`, is the natural
generalisation [Pancheshnyi2013]_ of the effective ionization coefficient
:math:`\alpha - \eta` to a chemistry with ion conversion and detachment.  It
is what :ref:`Chap:Eigenvalues` computes and tracks against :math:`E/N`.

This distinction matters.  Because :math:`\lambda_{\max}` accounts for
electrons that attach and are later released again, the field at which
:math:`\lambda_{\max} = 0` lies *below* the field at which
:math:`\alpha = \eta`.  A criterion built on :math:`\alpha - \eta` alone
therefore places the threshold too high in an electronegative gas, and the
gap between the two grows with :math:`pd` as the negative-ion transit time
becomes long enough for detachment to matter.
