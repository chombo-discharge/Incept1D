Theoretical model
=================

.. note::

   This page summarizes the derivation given in ``concepts.tex``, section
   "Methods" → "Theoretical model" (``\subsection{Theoretical model}``,
   ``\label{sec:model}``, starting around line 184). Equation labels below
   refer to the corresponding ``\label{...}`` in that file. Read the LaTeX
   source for the full derivation and literature references; this page only
   carries enough of it to map equations onto code.

Equations of motion
--------------------

The full three-dimensional evolution of electrons, ions, and ionizing
photons is approximated by a one-dimensional drift-reaction system. Let
:math:`\vec{n}(x, t)` be the column vector of species number densities
(length :math:`N_s`). The governing equation is (``eq:drift_reaction``)

.. math::

   \partial_t\vec{n} + \partial_x\left(\bm{V}\vec{n}\right)
     = \bm{R}\vec{n} + \sum_j\vec{\beta}_j\,\xi_j\,\kappa_j\,\Psi_j^0,

where :math:`\bm{V} = \mathrm{diag}(v_1, v_2, \ldots)` is the diagonal
matrix of signed species drift velocities and :math:`\bm{R}` is a
reaction-rate matrix that does **not** depend on :math:`\vec{n}` (local field
approximation: rates depend only on the reduced field :math:`E/N`,
temperature, and pressure). The last term is photoionization: :math:`\kappa_j
\Psi_j^0` is the number of ionizing photons of group :math:`j` absorbed per
unit volume and time, :math:`\xi_j` is a photoionization efficiency, and
:math:`\vec{\beta}_j` selects which species are produced by absorption of
that photon group.

Ionizing photons are assumed to be generated only by electrons, at a
volumetric rate :math:`\rho\,n_\mathrm{e}`. Rather than solving the full
radiative-transfer equation (which is only tractable in a slab geometry for
a genuinely 1-D model), a **two-stream approximation** is used: only photons
travelling along :math:`\pm x` are tracked, with a compensating
:math:`\Delta\Omega/(4\pi)` factor representing the fraction of generated
photons that are assumed to propagate along those two rays
(``eq:two_stream``):

.. math::

   \frac{1}{c}\partial_t\Psi_j^{\pm} \pm \partial_x\Psi_j^{\pm}
     = -\kappa_j\Psi_j^{\pm}
       + \left(\frac{\Delta\Omega}{4\pi}\right) g_j\,\rho\,\vec{u}_\mathrm{e}^\intercal \vec{n},

with :math:`\Psi_j^+(0) = 0` and :math:`\Psi_j^-(d) = 0` (no incoming
photons at either boundary), and :math:`\kappa_j\Psi_j^0 = \kappa_j(\Psi_j^+
+ \Psi_j^-)`. Charged-species diffusion is neglected throughout (bulk
avalanche growth, not sheath physics).

Code correspondence: the mechanism-file functions ``get_R``, ``get_V``,
``get_B``, ``get_C``, ``get_kappa`` supply :math:`\bm{R}`, :math:`\bm{V}`,
the photon-production coupling :math:`\bm{B} = [\vec{\beta}_1\xi_1\kappa_1,
\ldots]`, the photon-source coupling :math:`\bm{C} = \frac{\Delta\Omega}{4\pi}
\rho\,\vec{g}\,\vec{u}_\mathrm{e}^\intercal \bm{V}^{-1}`, and :math:`\bm{\kappa}
= \mathrm{diag}(\kappa_1, \kappa_2, \ldots)` respectively.

The augmented ODE and its formal solution
-------------------------------------------

Seeking exponential-in-time solutions :math:`\vec{n}(x,t) \to
\vec{n}(x)e^{\lambda t}`, :math:`\vec{\Psi}(x,t) \to \vec{\Psi}(x)e^{\lambda
t}` (valid because the coupled system is linear and homogeneous) turns the
PDE system into a linear ODE in :math:`x` for the augmented state vector
:math:`\vec{\theta} = (\vec{y}, \vec{\Psi}^+, \vec{\Psi}^-)^\intercal`, where
:math:`\vec{y} = \bm{V}\vec{n}` is the charged-species flux
(``eq:augmented_ode``):

.. math::

   \partial_x
   \begin{bmatrix} \vec{y} \\ \vec{\Psi}^+ \\ \vec{\Psi}^- \end{bmatrix}
   =
   \underbrace{
   \begin{bmatrix}
     \bm{A} & \bm{B} & \bm{B} \\
     \bm{C} & -\bm{D} & \bm{0} \\
     -\bm{C} & \bm{0} & \bm{D}
   \end{bmatrix}
   }_{\bm{\mathcal{A}}}
   \begin{bmatrix} \vec{y} \\ \vec{\Psi}^+ \\ \vec{\Psi}^- \end{bmatrix},
   \qquad
   \bm{A} = \left(\bm{R} - \lambda\bm{I}\right)\bm{V}^{-1},
   \qquad
   \bm{D} = \bm{\kappa} + \frac{\lambda}{c}\bm{I}.

This is exactly :func:`Inception._build_A_aug`: it assembles
:math:`\bm{\mathcal{A}}` (called ``A_aug`` in the code) for a given
:math:`E/N`, gap length, and :math:`\lambda`. (The code additionally
collapses photon groups whose optical depth across a sub-step is very large
into a *local* absorption term, to avoid propagating numerically stiff decay
modes — see the docstring of ``_build_A_aug`` for details; this is a
numerical convenience and does not change the physics.)

In a uniform gap :math:`\bm{\mathcal{A}}` is position-independent and the
formal solution is a matrix exponential, :math:`\vec{\theta}(x) =
e^{\bm{\mathcal{A}}x}\vec{\theta}_0 = \bm{M}(x)\vec{\theta}_0`
(``eq:theta_soln``). For a non-uniform field (sphere-plane / sphere-sphere
geometries, see :class:`FieldDistributions.FieldDistribution`)
:math:`\bm{\mathcal{A}}` varies with :math:`x` and :math:`\bm{M}(d)` is
instead built by composing per-segment propagators,
:math:`\bm{M}(d) = \prod_i \bm{P}_i`, where each :math:`\bm{P}_i` is a local
matrix exponential over a sub-interval (:func:`Inception.midpoint_propagator`,
:func:`Inception.magnus2_propagator`, with adaptive step refinement in
:func:`Inception._adaptive_midpoint_segment`).

Boundary conditions and the inception determinant
-----------------------------------------------------

Zero-flux conditions are imposed on negative ions at the cathode
(:math:`x=0`) and on positive ions at the anode (:math:`x=d`); photons have
open (no-influx) boundaries, :math:`\Psi_j^+(0)=0`, :math:`\Psi_j^-(d)=0`.
The cathode electron flux is a secondary-emission term sourced by ion and
photon bombardment (``eq:see_condition``):

.. math::

   \Pi_\mathrm{e}\vec{\theta}_0 = -\vec{\gamma}_\mathrm{p}^\intercal \Pi_\mathrm{p}\vec{\theta}_0
                                   + \vec{\gamma}_\Psi^\intercal \Pi_-\vec{\theta}_0,

where :math:`\Pi_\mathrm{e}`, :math:`\Pi_\mathrm{p}`, :math:`\Pi_\mathrm{m}`
select the electron, positive-ion, and negative-ion rows of
:math:`\vec{\theta}_0`, and :math:`\Pi_\pm` select the forward/backward
photon-flux rows. Collecting the cathode constraints (:math:`\bm{Q}_0
\vec{\theta}_0 = \vec{0}`) and the anode constraints, evaluated through the
propagator, :math:`\bm{Q}_d\vec{\theta}_0 = \vec{0}` (``eq:det_criterion``
and surrounding text) gives the full boundary-value problem

.. math::

   \bm{Q}(\lambda)\,\vec{\theta}_0 =
   \begin{bmatrix} \bm{Q}_0 \\ \bm{Q}_d \end{bmatrix} \vec{\theta}_0 = \vec{0}.

A non-trivial solution — i.e. a genuine discharge mode — exists if and only
if

.. math::

   \det \bm{Q}(\lambda) = 0.

:math:`\lambda > 0` is discharge growth, :math:`\lambda < 0` is decay, and
the **inception threshold** is :math:`\lambda = 0`. This determinant
condition is a field-line criterion for discharge inception that accounts
for primary and secondary ionization self-consistently, without requiring a
full 3-D discharge simulation.

Code correspondence:

- :func:`Inception._assemble_det_Q` builds :math:`\bm{Q}` from the
  propagator matrix :math:`\bm{M}(d)` and the mechanism's ``get_Pi_e`` /
  ``get_Pi_plus`` / ``get_Pi_minus`` / ``get_gamma_plus`` / ``get_gamma_Psi``,
  and returns :math:`\det\bm{Q}` (as a signed log-determinant with
  conditioning guards, returning ``NaN`` when :math:`\bm{Q}` is
  ill-conditioned).
- :func:`Inception.inception_det` is the end-to-end evaluation of
  :math:`\det\bm{Q}(\lambda)` for a given reference field, :math:`p\cdot d`,
  mechanism, and field geometry: it integrates :math:`\bm{\mathcal{A}}(x)`
  across the gap to build :math:`\bm{M}(d)`, then calls
  ``_assemble_det_Q``.
- :func:`Inception.find_all_breakdown_EN` brackets and root-finds
  :math:`E/N` such that ``inception_det(...) == 0`` at fixed :math:`p\cdot
  d`; :func:`Inception.compute_paschen_curve` repeats this over a
  :math:`p\cdot d` sweep and tracks multiple solution branches by
  continuity, producing a generalized Paschen curve.
- :func:`Lambda.find_lambda_for_voltage` instead fixes :math:`E/N` above
  threshold and solves for the temporal growth rate :math:`\lambda^* > 0`
  such that :math:`\det\bm{Q}(\lambda^*) = 0` — the rate at which the
  discharge current grows once the inception voltage is exceeded.

Recovering the standard (attachment-corrected) Paschen law
--------------------------------------------------------------

As a check, the general framework above reproduces the textbook result when
photoionization is dropped and the chemistry is reduced to three reactions:
impact ionization (:math:`\mathrm{e} \to 2\mathrm{e} + \mathrm{M}^+`,
coefficient :math:`\alpha`), attachment (:math:`\mathrm{e} \to
\mathrm{M}^-`, coefficient :math:`\eta`), and detachment (:math:`\mathrm{M}^-
\to \mathrm{M} + \mathrm{e}`, coefficient :math:`\delta`). The determinant
condition reduces to a closed form (``eq:generalized_paschen``):

.. math::

   \gamma\alpha\left[\frac{c_+}{\lambda_+}\left(e^{\lambda_+ d} - 1\right)
                      + \frac{c_-}{\lambda_-}\left(e^{\lambda_- d} - 1\right)\right] = 1,
   \qquad
   \lambda_\pm = \frac{\alpha - \eta - \delta \pm \Delta}{2},
   \qquad
   \Delta = \sqrt{(\alpha - \eta - \delta)^2 + 4\alpha\delta}.

:math:`\lambda_+` is the largest eigenvalue of :math:`\bm{R}\bm{V}^{-1}` and
is the *apparent effective ionization coefficient* — this is precisely what
:func:`Eigenvalues.compute_eigenvalues` computes and tracks vs. :math:`E/N`
for a full mechanism, generalized to more than three species. Setting
:math:`\delta = 0` recovers the classical attachment-corrected Paschen law

.. math::

   (\alpha - \eta) d = \ln\left(1 + \frac{1}{\gamma}\frac{\alpha - \eta}{\alpha}\right),

which further reduces to :math:`\alpha d = \ln(1 + \gamma^{-1})` when
:math:`\eta = 0`. Notably, solutions with :math:`\alpha \le \eta` also
exist (electrons attach faster than they ionize), interpreted via the
positive-ion/negative-ion growth rates rather than net electron growth; in
the limit :math:`\alpha \to \eta`, :math:`\alpha d = \gamma^{-1}`. This is
the regime referenced throughout the manuscript where "inception can occur
even when electron attachment exceeds electron impact ionization" once
detachment is restored (:math:`\delta \ne 0`).

A minimal reaction scheme for dry air
-----------------------------------------

The example mechanism (:mod:`Air.Air_Pancheshnyi`, ``Air/Air_Pancheshnyi.py``)
tracks electrons and five ion species
(:math:`\mathrm{N}_2^+, \mathrm{O}_2^+, \mathrm{O}^-, \mathrm{O}_2^-,
\mathrm{O}_3^-`) in 79%/21% N\ :sub:`2`/O\ :sub:`2` dry air, with reactions
grouped as:

- **Impact ionization** (reactions 1-2): BOLSIG+-computed rates feeding
  :math:`\alpha/N = (\nu_{\mathrm{N}_2} k_1 + \nu_{\mathrm{O}_2} k_2) /
  (\mu_\mathrm{e} E)`.
- **Attachment and detachment** (reactions 3-6): dissociative attachment to
  :math:`\mathrm{O}^-`, three-body attachment to :math:`\mathrm{O}_2^-`
  (Bloch-Bradbury mechanism, see below), collisional detachment from
  :math:`\mathrm{O}_2^-`, and associative detachment from :math:`\mathrm{O}^-`.
- **Ion conversion** (reactions 7-8): :math:`\mathrm{O}^- \to
  \mathrm{O}_2^-` and :math:`\mathrm{O}^- \to \mathrm{O}_3^-` (the latter is
  the only reaction in the scheme that acts as a *permanent* electron sink;
  :math:`\mathrm{O}_3^-` has a comparatively high electron-detachment energy
  barrier, :math:`\sim 2.1~\mathrm{eV}`).
- **Photoionization / photoemission** (reactions 9-10): the two-stream
  multigroup fit to the Zheleznyak absorption model
  (see ``Air/Zheleznyak.py`` and Appendix "Fitting the two-stream
  approximation to the Zheleznyak model" in ``concepts.tex``).
- **Secondary electron emission** (reactions 11-13): ion- and
  photon-induced SEE at the cathode.

The three-body attachment rate to :math:`\mathrm{O}_2^-` (reaction 4) is a
lumped representation of the sequential Bloch-Bradbury process (electron
capture → autodetachment / collisional stabilization); see the
``concepts.tex`` discussion around Eq. ``three_body_rate`` for when the
three-body approximation is (and is not) appropriate, which matters if you
extend the mechanism to pressures far outside the validated range
(the manuscript's regime is roughly :math:`p \in [10^{-3}, 10]~\mathrm{bar}`).

Reaction strings in the mechanism file are parsed declaratively by
:mod:`Reactions` (``"e + O2 + O2 -> O2- + O2"`` style strings); see
:doc:`modules/reactions` for the parsing rules.
