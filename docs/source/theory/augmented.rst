.. _Chap:AugmentedSystem:

The augmented system
====================

.. contents:: On this page
   :local:
   :depth: 1

Three ingredients are now in place: the charged-species ODE
:eq:`eq_flux_ode`, the two-stream photon equations :eq:`eq_two_stream`, and
the electrode conditions of :ref:`Chap:SecondaryEmission`.  They are coupled
— the photons are sourced by the electrons and in turn ionize the gas — so
they have to be solved together.  This page collects them into a single
first-order system and writes down its formal solution.

The augmented state vector
--------------------------

Collect the photon fluxes into vectors :math:`\vec{\Psi}^\pm = (\Psi_1^\pm,
\Psi_2^\pm, \ldots)^\intercal` and stack them under the charged-species
fluxes:

.. math::

   \vec{\theta} =
   \begin{bmatrix} \vec{y} \\ \vec{\Psi}^+ \\ \vec{\Psi}^- \end{bmatrix},
   \qquad \dim\vec{\theta} = N_s + 2N_\gamma .

This is the state of the gap at one position: what is flowing through
:math:`x`, and in which direction.

Block matrices
--------------

Define

.. math::
   :label: eq_blocks

   \begin{aligned}
   \bm{A} &= \left(\bm{R} - \lambda\bm{I}\right)\bm{V}^{-1}, &
   \dim\bm{A} &= N_s \times N_s, \\
   \bm{B} &= \begin{bmatrix}\vec{\beta}_1\xi_1\kappa_1 &
              \vec{\beta}_2\xi_2\kappa_2 & \cdots\end{bmatrix}, &
   \dim\bm{B} &= N_s \times N_\gamma, \\
   \bm{C} &= \frac{\Delta\Omega}{4\pi}\,\rho\,\vec{g}\,\vec{u}_\mathrm{e}^\intercal\bm{V}^{-1}, &
   \dim\bm{C} &= N_\gamma \times N_s, \\
   \bm{D} &= \bm{\kappa} + \frac{\lambda}{c}\bm{I}, \qquad
   \bm{\kappa} = \diag(\kappa_1, \kappa_2, \ldots), &
   \dim\bm{D} &= N_\gamma \times N_\gamma .
   \end{aligned}

Each block is one of the couplings already derived: :math:`\bm{A}` is the
plasma chemistry and drift of :eq:`eq_flux_ode`; :math:`\bm{B}` is photon
absorption producing charged species; :math:`\bm{C}` is electron impact
producing photons; and :math:`\bm{D}` is photon absorption plus the
retardation term :math:`\lambda/c`.

The augmented ODE
-----------------

With those, :eq:`eq_flux_ode` and :eq:`eq_two_stream` are one system:

.. math::
   :label: eq_augmented_ode

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
   \dim\bm{\mathcal{A}} = (N_s + 2N_\gamma)\times(N_s + 2N_\gamma).

The sign pattern in the lower blocks is the two-stream structure: the
forward and backward photon populations see the same source and the same
absorption, but propagate in opposite directions.

:math:`\bm{\mathcal{A}}` follows directly from the plasma chemistry and
photoionization data, and depends only on :math:`E/N`, :math:`p`, :math:`T`
and :math:`\lambda`.  When the mechanism has no photon groups
(:math:`N_\gamma = 0`) it reduces to :math:`\bm{A}` alone.

.. admonition:: Code

   :func:`incept1d.solver._build_A_aug` assembles :math:`\bm{\mathcal{A}}` for a
   given :math:`E/N`, :math:`p`, :math:`T`, and :math:`\lambda` from the
   mechanism's ``get_R``, ``get_V``, ``get_B``, ``get_C``, and ``get_kappa``.
   Note that the mechanism returns :math:`\bm{C}` *without* the
   :math:`\bm{V}^{-1}` factor (it acts on densities); the solver applies
   :math:`\bm{V}^{-1}` itself.  Every photon group is carried explicitly,
   however optically thick: replacing a thick group by the local term
   :math:`2\vec{b}_j\vec{c}_j^\intercal/\kappa_j` keeps the photoionization
   it produces but places each photoelectron at its point of emission, which
   removes the upstream seeding that makes photoionization a feedback loop.
   The stiffness of the :math:`e^{\pm\kappa_j x}` photon modes is handled by
   how the inception condition is evaluated
   (:ref:`Chap:Numerics:Riccati`).

Formal solution
---------------

The formal solution to :eq:`eq_augmented_ode` is

.. math::
   :label: eq_theta_soln

   \vec{\theta}(x) = \mathcal{P}\exp\left(\int_0^x\bm{\mathcal{A}}(s)\,ds\right)
   \vec{\theta}_0 = \bm{M}(x)\,\vec{\theta}_0,
   \qquad \vec{\theta}_0 = \vec{\theta}(0),

where :math:`\mathcal{P}\exp` is a path-ordered matrix exponential and
:math:`\bm{M}(x)` is the **propagator** from the cathode to :math:`x`.

* In a **uniform gap** :math:`\bm{\mathcal{A}}` is position-independent and
  :math:`\bm{M}(x) = e^{\bm{\mathcal{A}}x}` is an ordinary matrix
  exponential.
* In a **non-uniform gap** :math:`\bm{\mathcal{A}}` varies with :math:`x`
  through :math:`E(x)`.  The domain is decomposed into sub-intervals on an
  adaptive grid, a local propagator :math:`\bm{P}_i` is computed on each
  (midpoint rule or a two-term Magnus expansion with Gauss-Legendre
  quadrature), and the total propagator is the path-ordered product
  :math:`\bm{M}(d) = \bm{P}_N\cdots\bm{P}_2\bm{P}_1`.  Details are in
  :ref:`Chap:Numerics:Propagator`.

Everything that follows only requires :math:`\bm{M}(d)`.  The anode
conditions at :math:`x = d` are pulled back to :math:`x = 0` through it, and
combined with the cathode block :math:`\bm{Q}_0` of :eq:`eq_Q0` the
inception problem becomes a single linear system for :math:`\vec{\theta}_0`
alone — which is :ref:`Chap:InceptionCriterion`.
