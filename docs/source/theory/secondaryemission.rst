.. _Chap:SecondaryEmission:

Secondary emission and approximations
=====================================

.. contents:: On this page
   :local:
   :depth: 1

Boundary conditions at the electrodes
-------------------------------------

The augmented ODE :eq:`eq_augmented_ode` is a first-order system in
:math:`x` on the domain :math:`[0, d]`, so :math:`N_s + 2N_\gamma`
conditions are needed to pin down :math:`\vec{\theta}_0 = \vec{\theta}(0)`.
They are split between the two electrodes:

At the **cathode** (:math:`x = 0`):

* Zero incoming flux of negative ions — they drift *away* from the cathode,
  so nothing enters;
* Zero incoming forward photon flux, :math:`\Psi_j^+(0) = 0`;
* The electron flux is the **secondary emission** flux produced by ions and
  photons hitting the cathode (below).

At the **anode** (:math:`x = d`):

* Zero incoming flux of positive ions, :math:`\Pi_+\vec{\theta}(d) = \vec{0}`;
* Zero incoming backward photon flux, :math:`\Psi_j^-(d) = 0`.

Counting: :math:`N_-` negative-ion conditions, :math:`1` electron condition,
and :math:`N_\gamma` photon conditions at the cathode; :math:`N_+`
positive-ion conditions and :math:`N_\gamma` photon conditions at the anode.
With :math:`N_s = 1 + N_+ + N_-` this is exactly :math:`N_s + 2N_\gamma`
conditions.  Every tracked species must therefore be classified as the
electron, a positive ion, or a negative ion — this is what the row-selection
matrices below encode.

Row-selection matrices
----------------------

Let :math:`\Pi_\mathrm{e}`, :math:`\Pi_+`, and :math:`\Pi_-` be the
row-selection matrices that pick the electron row, the positive-ion rows,
and the negative-ion rows of :math:`\vec{\theta}` (:math:`\Pi_\mathrm{e}` is
a single row).  Similarly let :math:`\Pi_\rightarrow` and
:math:`\Pi_\leftarrow` select the :math:`\vec{\Psi}^+` and
:math:`\vec{\Psi}^-` blocks.  For the six-species dry-air scheme
(:math:`\mathrm{e}, \mathrm{N}_2^+, \mathrm{O}_2^+, \mathrm{O}^-,
\mathrm{O}_2^-, \mathrm{O}_3^-`):

.. math::

   \Pi_\mathrm{e} = \begin{bmatrix}1&0&0&0&0&0\end{bmatrix},\quad
   \Pi_+ = \begin{bmatrix}0&1&0&0&0&0\\0&0&1&0&0&0\end{bmatrix},\quad
   \Pi_- = \begin{bmatrix}0&0&0&1&0&0\\0&0&0&0&1&0\\0&0&0&0&0&1\end{bmatrix}.

.. admonition:: Code

   A mechanism file returns these through ``get_Pi_e()``, ``get_Pi_plus()``
   and ``get_Pi_minus()``, each of shape ``(rows, N_s)``.  The photon
   selectors are built by the solver, since it decides at run time which
   photon groups are propagated explicitly (:ref:`Chap:Numerics:Propagator`).

The secondary-emission condition
--------------------------------

Let :math:`\vec{\gamma}_\mathrm{p} = \vec{\gamma}_\mathrm{p}(E)` be the
column vector of SEE efficiencies due to ion bombardment (one entry per
positive-ion species, in the row order of :math:`\Pi_+`), and
:math:`\vec{\gamma}_\Psi` the efficiencies due to photon bombardment (one
entry per photon group).  The emitted electron flux is

.. math::
   :label: eq_see_condition

   \Pi_\mathrm{e}\vec{\theta}_0 =
     -\vec{\gamma}_\mathrm{p}^\intercal\,\Pi_+\vec{\theta}_0
     + \vec{\gamma}_\Psi^\intercal\,\Pi_\leftarrow\vec{\theta}_0 .

The minus sign on the ion term comes from the coordinate convention: the
positive-ion flux at the cathode is negative (ions move towards
:math:`-x`), while the emitted electron flux is positive.  The backward
photon flux :math:`\vec{\Psi}^-(0)` is positive by definition.

Combining the three cathode constraints into a block system
:math:`\bm{Q}_0\vec{\theta}_0 = \vec{0}` gives

.. math::
   :label: eq_Q0

   \bm{Q}_0 = \begin{bmatrix}
     \Pi_\mathrm{e} + \vec{\gamma}_\mathrm{p}^\intercal\Pi_+
                    - \vec{\gamma}_\Psi^\intercal\Pi_\leftarrow \\
     \Pi_- \\
     \Pi_\rightarrow
   \end{bmatrix} .

Ion-induced emission
--------------------

The ion-induced yields are the principal free parameters of the cathode
model.  The dry-air mechanism uses a field-dependent form

.. math::

   \gamma_i(E) = \gamma_0 + \gamma_1\exp\left(-\frac{E_\mathrm{ref}}{\beta E}\right),
   \qquad E = (E/N)\,N ,

identical for :math:`\mathrm{N}_2^+` and :math:`\mathrm{O}_2^+`, with the
default :math:`\gamma_0 = 10^{-3}`, :math:`\gamma_1 = 0` (i.e.
a constant yield).  The four parameters are exposed as ``gamma0``,
``gamma1``, ``eref``, ``beta`` in the JSON configuration and can be set
independently for the two polarities of a non-symmetric gap
(:ref:`Chap:Configuration`).  There is considerable uncertainty in these
constants; ``mechanisms/air/see.json`` sweeps :math:`\gamma_0` over
:math:`10^{-4}`–:math:`10^{-2}` to quantify the sensitivity.

.. admonition:: Code

   ``get_gamma_plus(EN, p, T)`` returns :math:`\vec{\gamma}_\mathrm{p}` with
   the module defaults; ``get_gamma_plus_with(EN, p, T, gamma0=, gamma1=,
   eref=, beta=)`` is the same with per-call overrides.  The solver always
   calls the latter (through :class:`incept1d.mechanism.Mechanism`), which is how
   polarity-specific overrides are applied without reloading the mechanism.
   The yields are evaluated at the *cathode* field, which for a non-uniform
   gap differs between polarities.

Photon-induced emission
-----------------------

Photoemission by ionizing photons is treated with a constant yield
:math:`\gamma_\Psi = 0.1` for every group in the dry-air scheme.  Because
the model only tracks photons energetic enough to ionize O\ :sub:`2`
(:math:`> 12` eV) and not the broader spectrum down to the cathode work
function, this term is likely an underestimate of the true photoemission
feedback.  The overall strength
of the photon feedback can be scaled with ``xi_emit`` (photoemission) and
``xi_photo`` (photoionization) in the JSON configuration, and switched off
with ``cone_angle: 0``.

.. admonition:: Code

   ``get_gamma_Psi(EN, p, T)`` returns :math:`\vec{\gamma}_\Psi`, shape
   :math:`(N_\gamma,)`.  :class:`incept1d.mechanism.Mechanism` multiplies it by
   ``xi_emit`` and the output of ``get_B`` by ``xi_photo``.
