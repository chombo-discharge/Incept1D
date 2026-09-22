.. _Chap:Transport:

Transport equations and approximations
======================================

.. contents:: On this page
   :local:
   :depth: 1

Coordinate convention
---------------------

The spatial coordinate :math:`x` runs from the **cathode at** :math:`x = 0`
**to the anode at** :math:`x = d`.  The electric field points in the
:math:`-x` direction, but the code works exclusively with the magnitude
:math:`E/N > 0`.  A mechanism file encodes the sign of each species' motion
in the drift-velocity matrix,

.. math::

   V_{ii} = -\operatorname{sgn}(Z_i)\,\mu_i\,|E| ,

so that electrons and negative ions (:math:`Z_i = -1`) have positive
velocity (towards the anode) and positive ions (:math:`Z_i = +1`) have
negative velocity (towards the cathode).  With this convention a positive
eigenvalue of :math:`\bm{R}\bm{V}^{-1}` corresponds to a spatially growing
solution — see :ref:`Chap:Eigenvalues`.

For non-uniform gaps the geometry is described by a normalised profile
:math:`f(\xi)`, :math:`\xi = x/d \in [0,1]`, with :math:`\int_0^1 f\,d\xi = 1`,
so that :math:`E(x) = E_\mathrm{ref}\,f(x/d)` and :math:`E_\mathrm{ref} =
U/d` is the mean field (:ref:`Chap:FieldDistributions`).  The profile is
defined with :math:`\xi = 0` at the high-field electrode; for a sphere-plane
gap the solver evaluates both polarities by flipping the profile.

The drift-reaction system
-------------------------

Starting from :eq:`eq_drift_reaction`, follow the eigenvalue approach of
[Ferreira2019]_ and seek solutions with exponential time dependence, :math:`\vec{n}(x,t) \to \vec{n}(x)e^{\lambda t}` and
:math:`\vec{\Psi}(x,t) \to \vec{\Psi}(x)e^{\lambda t}`.  This is allowed
without loss of generality because the system is linear and homogeneous.
Writing the equations in terms of the **charged-species flux**
:math:`\vec{y} = \bm{V}\vec{n}` rather than the density gives

.. math::
   :label: eq_flux_ode

   \partial_x\vec{y} = \left(\bm{R} - \lambda\bm{I}\right)\bm{V}^{-1}\vec{y}
     + \sum_j\vec{\beta}_j\xi_j\kappa_j\left(\Psi_j^+ + \Psi_j^-\right),

where :math:`\Psi_j^\pm` are the forward/backward photon fluxes of group
:math:`j` introduced in :ref:`Chap:Photoionization`.  Working with fluxes
is what makes the boundary conditions simple: zero-flux conditions on ions
and secondary-emission conditions on the cathode are linear constraints on
:math:`\vec{y}` at a point.

The augmented ODE
-----------------

Collect the photon fluxes into vectors :math:`\vec{\Psi}^\pm = (\Psi_1^\pm,
\Psi_2^\pm, \ldots)^\intercal` and define the block matrices

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

Then :eq:`eq_flux_ode` and the two-stream equations
:eq:`eq_two_stream` combine into a single first-order linear ODE for the
augmented state vector :math:`\vec{\theta} = (\vec{y}, \vec{\Psi}^+,
\vec{\Psi}^-)^\intercal`:

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

:math:`\bm{\mathcal{A}}` is trivially definable from the plasma chemistry
and photoionization data and depends only on :math:`E/N`, :math:`p`,
:math:`T`, and :math:`\lambda`.  When the mechanism has no photon groups
(:math:`N_\gamma = 0`) it reduces to :math:`\bm{A}` alone.

.. admonition:: Code

   :func:`incept1d.solver._build_A_aug` assembles :math:`\bm{\mathcal{A}}` for a
   given :math:`E/N`, :math:`p`, :math:`T`, and :math:`\lambda` from the
   mechanism's ``get_R``, ``get_V``, ``get_B``, ``get_C``, and ``get_kappa``.
   Note that the mechanism returns :math:`\bm{C}` *without* the
   :math:`\bm{V}^{-1}` factor (it acts on densities); the solver applies
   :math:`\bm{V}^{-1}` itself.  Photon groups whose optical depth
   :math:`\kappa_j d` across the gap exceeds a threshold are not propagated
   explicitly but folded into :math:`\bm{A}` as a local absorption term
   :math:`2\vec{b}_j\vec{c}_j^\intercal/\kappa_j`; see :ref:`Chap:Numerics:Propagator`.

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

Everything that follows only requires :math:`\bm{M}(d)`: the boundary
conditions at :math:`x = d` are pulled back to :math:`x = 0` through it,
and the inception problem becomes a linear system for :math:`\vec{\theta}_0`
alone (:ref:`Chap:InceptionCriterion`).

The transport matrix and its eigenvalues
----------------------------------------

At :math:`\lambda = 0` and without photons, :eq:`eq_augmented_ode` reduces to
:math:`\partial_x\vec{y} = \bm{R}\bm{V}^{-1}\vec{y}`.  The eigenvalues
:math:`k_i` of :math:`\bm{R}\bm{V}^{-1}` are the *local* spatial growth
rates of the flux modes: the solution can be written as a sum of eigenmodes
:math:`\sum_i c_i\vec{l}_ie^{k_ix}`, and if any :math:`\mathrm{Re}\,k_i > 0`
there is net spatial growth.  The largest such eigenvalue,
:math:`\lambda_{\max}`, is the natural generalisation [Pancheshnyi2013]_ of the effective
ionization coefficient :math:`\alpha - \eta` to a chemistry with ion
conversion and detachment, and is what :ref:`Chap:Eigenvalues` computes and
tracks as a function of :math:`E/N`.  Because it includes detachment, the
field at which :math:`\lambda_{\max} = 0` lies *below* the field at which
:math:`\alpha = \eta` — this is the origin of the reduced breakdown voltages
at large :math:`pd` (see :ref:`Chap:AirScheme`).
