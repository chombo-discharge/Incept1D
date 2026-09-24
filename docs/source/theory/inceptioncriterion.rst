.. _Chap:InceptionCriterion:

Inception criterion
===================

.. contents:: On this page
   :local:
   :depth: 1

The boundary-value problem
--------------------------

The cathode conditions were collected into :math:`\bm{Q}_0\vec{\theta}_0 =
\vec{0}` in :eq:`eq_Q0`.  The anode conditions are imposed on
:math:`\vec{\theta}(d)`, which the formal solution :eq:`eq_theta_soln`
expresses through the propagator as :math:`\vec{\theta}(d) =
\bm{M}(d)\vec{\theta}_0`:

.. math::
   :label: eq_Qd

   \bm{\Pi}_+\bm{M}(d)\vec{\theta}_0 = \vec{0}, \qquad
   \bm{\Pi}_{\Psi^-}\bm{M}(d)\vec{\theta}_0 = \vec{0}
   \quad\Longrightarrow\quad
   \bm{Q}_d = \begin{bmatrix}\bm{\Pi}_+\bm{M}(d)\\ \bm{\Pi}_{\Psi^-}\bm{M}(d)\end{bmatrix} .

The complete specification of the boundary-value problem is the square,
homogeneous linear system

.. math::
   :label: eq_Q_system

   \bm{Q}(\lambda)\,\vec{\theta}_0 =
   \begin{bmatrix}\bm{Q}_0\\ \bm{Q}_d\end{bmatrix}\vec{\theta}_0 = \vec{0},
   \qquad \dim\bm{Q} = (N_s + 2N_\gamma)\times(N_s + 2N_\gamma).

The determinant criterion
-------------------------

A non-trivial solution of :eq:`eq_Q_system` — a self-consistent discharge
mode — exists if and only if

.. math::
   :label: eq_det_criterion

   \det\bm{Q}(\lambda) = 0 .

Since :math:`\bm{Q}` depends on the applied field through
:math:`\bm{M}(d)` and :math:`\vec{\gamma}_+`, and on :math:`\lambda`
through :math:`\bm{A}` and :math:`\bm{D}`, :eq:`eq_det_criterion` is a
relation between the applied voltage and the temporal growth rate:

* :math:`\lambda > 0` indicates discharge **growth**, :math:`\lambda < 0`
  **decay**.
* The **inception threshold** is :math:`\lambda = 0`: the voltage
  :math:`U^*` (equivalently :math:`(E/N)^*` at fixed :math:`pd`) at which
  :math:`\det\bm{Q}(0) = 0`.  This is what :ref:`Chap:Inception` solves for
  along a :math:`pd` sweep, producing an inception curve (a generalized Paschen curve).
* Fixing the voltage above :math:`U^*` and treating :math:`\lambda` as the
  unknown gives the **temporal growth rate** of the discharge, which is what
  :ref:`Chap:Lambda` computes.

:eq:`eq_det_criterion` is, in essence, a field-line criterion for the
inception of electrical discharges that accounts for primary and secondary
ionization self-consistently.  Compared to fully three-dimensional discharge
simulations it is computationally trivial: one determinant per
:math:`(E/N, pd)` evaluation, with :math:`\bm{M}(d)` obtained from a handful
of small matrix exponentials.

.. admonition:: Code

   :func:`incept1d.solver._det_Q` builds :math:`\bm{Q}` from
   :math:`\bm{M}(d)`, the mechanism's row selectors and SEE yields, and
   returns :math:`\det\bm{Q}` (row-normalised, via a signed
   log-determinant, with ``NaN`` returned when :math:`\bm{Q}` is
   ill-conditioned, or evaluated without forming :math:`\bm{M}` if
   requested).  :func:`incept1d.solver.inception_det` is the end-to-end
   evaluation for given :math:`E/N`, :math:`pd`, :math:`\lambda`, and field
   geometry.  Root-finding in :math:`E/N` is done by
   :func:`incept1d.inception.find_all_breakdown_EN`; root-finding in :math:`\lambda`
   by :func:`Lambda.find_lambda_for_voltage`.  See :ref:`Chap:Numerics:RootFinding`
   for how sign changes and ``NaN`` regions are handled.

Recovering the standard Paschen law
-----------------------------------

The general framework reproduces the textbook result.  Let :math:`\mathrm{M}`
denote a neutral molecule, drop photoionization, and describe the charged
species by three reactions: impact ionization
:math:`\mathrm{e} \xrightarrow{k_1} 2\mathrm{e} + \mathrm{M}^+`, attachment
:math:`\mathrm{e} \xrightarrow{k_2} \mathrm{M}^-`, and detachment
:math:`\mathrm{M}^- \xrightarrow{k_3} \mathrm{M} + \mathrm{e}`.  With
:math:`\vec{y} = (y_\mathrm{e}, y_+, y_-)` and the Townsend coefficients
:math:`\alpha = k_1/v_\mathrm{e}`, :math:`\eta = k_2/v_\mathrm{e}`,
:math:`\delta = k_3/v_-`, the inception matrix at :math:`\lambda = 0` is

.. math::

   \bm{\mathcal{A}} = \bm{R}\bm{V}^{-1} =
   \begin{bmatrix}
     \alpha - \eta & 0 & \delta \\
     \alpha & 0 & 0 \\
     \eta & 0 & -\delta
   \end{bmatrix},

with :math:`\bm{\Pi}_\mathrm{e} = (1,0,0)`, :math:`\bm{\Pi}_+ = (0,1,0)`,
:math:`\bm{\Pi}_- = (0,0,1)` and a scalar yield :math:`\gamma`.  The matrix
exponential can be computed analytically, and :math:`\det\bm{Q}(0) = 0`
becomes

.. math::
   :label: eq_generalized_paschen

   \gamma\alpha\left[\frac{c_+}{\lambda_+}\left(e^{\lambda_+ d} - 1\right)
                   + \frac{c_-}{\lambda_-}\left(e^{\lambda_- d} - 1\right)\right] = 1,

with

.. math::

   \lambda_\pm = \frac{\alpha - \eta - \delta \pm \Delta}{2}, \qquad
   c_\pm = \frac{\Delta \pm (\alpha - \eta + \delta)}{2\Delta}, \qquad
   \Delta = \sqrt{(\alpha - \eta - \delta)^2 + 4\alpha\delta}.

Here :math:`\lambda_\pm` are the two eigenvalues of the reduced model's
electron/negative-ion block, so :math:`\lambda_+` is its
:math:`\lambda_{\max}` in the sense of :ref:`Chap:Transport`: the largest
eigenvalue of
:math:`\bm{R}\bm{V}^{-1}` and defines the *apparent effective ionization
coefficient*; it reduces to :math:`\alpha - \eta` when :math:`\delta = 0`.
It is :math:`\lambda_+`, not :math:`\alpha - \eta`, that determines the net
growth of charge carriers through the gap, so the critical field at which
growth ceases is set by the detachment kinetics and not by
:math:`\alpha = \eta` alone.

When :math:`\delta = 0`, :eq:`eq_generalized_paschen` reduces to the
attachment-corrected Paschen law

.. math::
   :label: eq_standard_paschen

   \gamma\,\frac{\alpha}{\alpha - \eta}\left(e^{(\alpha - \eta)d} - 1\right) = 1
   \quad\Longleftrightarrow\quad
   (\alpha - \eta)\,d = \ln\left(1 + \frac{1}{\gamma}\frac{\alpha - \eta}{\alpha}\right),

which for :math:`\eta = 0` or :math:`\alpha \gg \eta` becomes the familiar
:math:`\alpha d = \ln(1 + \gamma^{-1})`.

Solutions in the attachment regime
..................................

:math:`\alpha > \eta` is *not* required for a solution.
:eq:`eq_standard_paschen` also has solutions when :math:`\alpha \le \eta`,
which requires a slight reinterpretation of the main charge carriers.  When
:math:`\eta > \alpha`, electrons released from the cathode are more likely
to attach than to ionize and the electron flux decays as
:math:`e^{(\alpha - \eta)x}`, but the cathode-emitted electrons still
produce positive ions before their avalanches terminate.  For
:math:`(\eta - \alpha)d \gg 1` the number of positive ions produced per
emitted electron is :math:`\alpha/(\eta - \alpha)`, and the number of SEE
electrons is :math:`\gamma\alpha/(\eta - \alpha)`, which may exceed one.
In the limit :math:`\alpha \to \eta` the criterion becomes
:math:`\alpha d = \gamma^{-1}`: although there is no net growth of the
electron flux, positive ions still grow at rate :math:`\alpha` and negative
ions at rate :math:`\eta`.  With common values of :math:`\gamma` this means
the gap must exceed a few hundred to a few thousand avalanche lengths.

Beyond these elementary considerations the textbook forms remain limited:
:eq:`eq_generalized_paschen` does not account for conversion between
different negative-ion species, only a subset of which may be unstable.
In air, dissociative attachment produces :math:`\mathrm{O}^-`, which can
later convert to :math:`\mathrm{O}_2^-` (unstable) or :math:`\mathrm{O}_3^-`
(comparatively stable) — precisely the situation the full determinant
criterion handles and the closed form does not.

.. admonition:: Code

   The :math:`3\times3` reduction is the analytic reference case for
   validating changes to the propagator or determinant code: a mechanism
   with only ionization, attachment, and detachment must reproduce
   :eq:`eq_generalized_paschen`, and with detachment switched off it must
   reproduce :eq:`eq_standard_paschen`.  ``mechanisms/air/pancheshnyi/paschen.json`` configures
   the dry-air mechanism into this limit (no detachment, no ion conversion,
   no photon feedback).

Derived criteria: ionization integrals
--------------------------------------

To interpret the inception mechanism it is useful to compare the full
criterion against two purely local surrogates evaluated along the field
line.  The **effective ionization integral** is

.. math::
   :label: eq_ionization_integral

   I_\alpha = \int_0^d \max\left(\alpha - \eta, 0\right)dx,

and the **apparent effective ionization integral** is

.. math::
   :label: eq_apparent_ionization_integral

   I_\lambda = \int_0^d \max\left(\lambda_{\max}, 0\right)dx,

where :math:`\lambda_{\max}` is the largest eigenvalue of
:math:`\bm{R}\bm{V}^{-1}`.  :math:`I_\alpha` counts only impact ionization
as an electron source, whereas :math:`I_\lambda` also accounts for ion
conversion, transport, detachment, and photoionization.  The **streamer
criterion** :math:`I_\alpha = C` with :math:`C \approx 18` is the customary
engineering inception criterion for non-uniform fields; the code can solve
for the voltage that satisfies it and overlay it on the inception curve.

.. admonition:: Code

   :func:`IonizationIntegral.aed_integral` and
   :func:`IonizationIntegral.eig_integral` evaluate :math:`I_\alpha` and
   :math:`I_\lambda`; ``incept1d pdiv --streamer-criterion C`` solves
   :math:`I_\alpha = C`.  Both require the mechanism to expose ``alpha`` and
   ``eta`` functions in addition to the required interface.
