.. _Chap:Photoionization:

Photoionization and approximations
==================================

.. contents:: On this page
   :local:
   :depth: 1

Radiative transfer
------------------

Ionizing photons are assumed to be generated only by electron impact on
neutral molecules, at a volumetric rate :math:`\rho\,n_\mathrm{e}` where
:math:`n_\mathrm{e} = \vec{u}_\mathrm{e}^\intercal\vec{n}` and
:math:`\vec{u}_\mathrm{e}` is the electron unit vector in species space.
Assuming isotropic emission, the radiative transfer equation for photon
group :math:`j` is

.. math::

   \frac{1}{c}\partial_t\Psi_j + \hat{\bm{\Omega}}\cdot\nabla\Psi_j
     = -\kappa_j\Psi_j + \frac{1}{4\pi}g_j\,\rho\,\vec{u}_\mathrm{e}^\intercal\vec{n},

where :math:`\Psi_j` is the photon flux in direction
:math:`\hat{\bm{\Omega}}`, :math:`\kappa_j` is the absorption coefficient of
the group, and :math:`g_j` (with :math:`\sum_j g_j = 1`) is the fraction of
photons emitted into group :math:`j`.  Rather than a Helmholtz/Eddington
closure [Bourdon2007]_ — for which consistent boundary conditions are
difficult to impose [Marskar2019]_ — the code works directly from a
truncated discrete-ordinates form of this equation.

The two-stream approximation
----------------------------

A consistent one-dimensional model of the genuinely three-dimensional
radiative transfer is only possible in a slab geometry (infinite transverse
extent), which is not realistic for a developing discharge.  To avoid
solving a full 3-D problem, a **two-stream approximation** is imposed: only
photons travelling along :math:`\pm x` are tracked, and photons propagating
in other directions are ignored.

In a uniform field, transverse photons would mostly seed electrons beside
the developing discharge and lead to broadening; in a strongly non-uniform
field they may cause photoionization on neighbouring field lines that later
curve back into the inception region.  Neither effect is captured.  What
*is* captured is the amount of photon feedback along the field line,
including photoemission at the cathode.

Assuming that *all* generated photons travel along the two rays would
overestimate the ionizing flux, so a conical angle :math:`\Delta\Omega` is
introduced: only photons emitted into this solid angle are assumed to
contribute.  The two-stream equations with the compensated source are

.. math::
   :label: eq_two_stream

   \frac{1}{c}\partial_t\Psi_j^{\pm} \pm \partial_x\Psi_j^{\pm}
     = -\kappa_j\Psi_j^{\pm}
       + \left(\frac{\Delta\Omega}{4\pi}\right) g_j\,\rho\,\vec{u}_\mathrm{e}^\intercal\vec{n},

with :math:`\Psi_j^+` the flux in :math:`+x` (towards the anode) and
:math:`\Psi_j^-` the flux in :math:`-x` (towards the cathode).  The number
of photons absorbed per unit volume and time — the coupling term in
:eq:`eq_drift_reaction` — is :math:`\kappa_j\Psi_j^0 = \kappa_j(\Psi_j^+ +
\Psi_j^-)`.  The boundary conditions are open, with no incoming photons:

.. math::

   \Psi_j^+(0) = 0, \qquad \Psi_j^-(d) = 0 .

Cone angle
..........

The fraction :math:`\Delta\Omega/(4\pi)` is a free parameter with
substantial geometric uncertainty.  In the code it is set from a
half-opening angle :math:`\theta_\mathrm{cone}` via

.. math::

   \frac{\Delta\Omega}{4\pi} = \frac{1 - \cos\theta_\mathrm{cone}}{2},

so :math:`\theta_\mathrm{cone} = 45^\circ` (the default) gives
:math:`\Delta\Omega/(4\pi) \approx 0.146`, :math:`90^\circ` gives
:math:`1/2` (a full hemisphere), and :math:`0^\circ` switches photon
feedback off entirely.  It is exposed as ``cone_angle`` in the JSON
configuration (:ref:`Chap:Configuration`); ``mechanisms/air/see.json`` contains a
ready-made sensitivity sweep over it.

After the eigenvalue substitution :math:`\vec{\Psi}(x,t) \to
\vec{\Psi}(x)e^{\lambda t}`, :eq:`eq_two_stream` becomes

.. math::

   \partial_x\vec{\Psi}^\pm = \mp\left(\bm{\kappa} + \frac{\lambda}{c}\bm{I}\right)\vec{\Psi}^\pm
     \pm \frac{\Delta\Omega}{4\pi}\,\rho\,\vec{g}\,\vec{u}_\mathrm{e}^\intercal\bm{V}^{-1}\vec{y},

which supplies the :math:`\bm{C}` and :math:`\bm{D}` blocks of the
augmented matrix :eq:`eq_augmented_ode`, while the absorption term
:math:`\sum_j\vec{\beta}_j\xi_j\kappa_j(\Psi_j^+ + \Psi_j^-)` in
:eq:`eq_flux_ode` supplies :math:`\bm{B}`.  The :math:`\lambda/c` term is
retardation; it is negligible at inception (:math:`\lambda = 0`) but is kept
for the growth-rate calculation in :ref:`Chap:Lambda`.

.. admonition:: Code

   A mechanism file supplies the photon data through four functions
   (see :ref:`Chap:NewMechanisms`):

   * ``get_kappa(p, T)`` → :math:`(\kappa_1, \kappa_2, \ldots)`, shape
     :math:`(N_\gamma,)`, in :math:`\mathrm{m}^{-1}`;
   * ``get_B(EN, p, T)`` → :math:`\bm{B}`, shape :math:`(N_s, N_\gamma)`,
     with :math:`B_{ij} = \beta_j[i]\,\xi_j\,\kappa_j`;
   * ``get_C(EN, p, T)`` → :math:`\frac{\Delta\Omega}{4\pi}\rho\,\vec{g}\,\vec{u}_\mathrm{e}^\intercal`,
     shape :math:`(N_\gamma, N_s)`, acting on *densities*; the solver
     multiplies by :math:`\bm{V}^{-1}` (:func:`incept1d.solver._build_A_aug`);
   * ``get_gamma_Psi(EN, p, T)`` → photoemission yields, see
     :ref:`Chap:SecondaryEmission`.

   A mechanism without photoionization returns zero-width arrays
   (:math:`N_\gamma = 0`).

Where the numbers come from
---------------------------

Nothing above fixes :math:`\rho`, :math:`g_j` or :math:`\kappa_j`: the
two-stream model is a transport scheme, and the spectrum it transports is a
property of the gas.  A mechanism must supply them, by whatever route suits
it — a measured absorption curve, a fitted multigroup decomposition, or a
single effective group.

The dry-air mechanism shipped with the code fits its groups to the
Zheleznyak absorption function; that construction, and the numbers it
produces, are in :ref:`Chap:AirScheme`.
