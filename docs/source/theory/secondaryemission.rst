.. _Chap:SecondaryEmission:

Secondary emission and approximations
=====================================

.. contents:: On this page
   :local:
   :depth: 1

Boundary conditions at the electrodes
-------------------------------------

The gap carries :math:`N_s` charged-species fluxes :math:`\vec{y}`
(:ref:`Chap:Transport`) and :math:`2N_\gamma` photon fluxes
:math:`\vec{\Psi}^\pm` (:ref:`Chap:Photoionization`), each governed by a
first-order equation in :math:`x`.  Pinning the solution down therefore takes
:math:`N_s + 2N_\gamma` conditions, split between the two electrodes.

Each one says the same kind of thing: what crosses that electrode, and in
which direction.

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

A positive ion arriving at the cathode may liberate an electron from the
surface.  The probability that it does is
:math:`\vec{\gamma}_\mathrm{p} = \vec{\gamma}_\mathrm{p}(E)` — one number
per positive-ion species, in general a function of the field at the cathode,
since a faster ion arrives with more energy.

That probability is an **input to the model, not a result of it**.  The
solver asks the mechanism for it and uses whatever comes back; it makes no
assumption about the functional form, and none is needed for the derivation.
A mechanism is free to return a constant, a fitted field dependence, or a
table.

This is worth stating plainly because :math:`\vec{\gamma}_\mathrm{p}` is
the least certain quantity in the whole model.  It depends on the cathode
material, its oxide layer, its roughness and its history, and quoted values
for the same nominal surface span orders of magnitude.  The criterion is only
logarithmically sensitive to it — :math:`\ln(1 + \gamma^{-1})` in the
classical limit — which is what makes the approach usable at all, but a
sensitivity sweep is still the honest way to report a result.  See
:ref:`Chap:AirScheme` for the form and values the dry-air mechanism uses, and
``mechanisms/air/see.json`` for a ready-made sweep.

.. admonition:: Code

   ``get_gamma_plus(EN, p, T)`` returns :math:`\vec{\gamma}_\mathrm{p}` with
   the mechanism's own defaults; ``get_gamma_plus_with(EN, p, T, ...)`` is the
   same with per-call overrides.  The solver always calls the latter (through
   :class:`incept1d.mechanism.Mechanism`), which is how polarity-specific
   overrides are applied without reloading the mechanism.  The yields are
   evaluated at the *cathode* field, which for a non-uniform gap differs
   between polarities.

The photoelectric effect
------------------------

The second way the cathode returns an electron is optical.  An ionizing
photon produced in the gas can travel *back* to the cathode and eject an
electron from the surface — the photoelectric effect — with probability
:math:`\vec{\gamma}_\Psi`, one entry per photon group.  This is the
:math:`\vec{\gamma}_\Psi^\intercal\Pi_\leftarrow` term in
:eq:`eq_see_condition`, acting on the backward flux
:math:`\vec{\Psi}^-(0)`.

It matters for a reason worth being explicit about: it is *fast*.  Ion-induced
emission requires an ion to drift the length of the gap, which at atmospheric
pressure takes microseconds.  A photon crosses the same gap in nanoseconds.
Where photoemission is strong enough to sustain the discharge on its own, the
cathode feedback loop closes many orders of magnitude sooner, and the relevant
question stops being *whether* the gap breaks down and becomes *how fast* —
which is what :ref:`Chap:Lambda` computes.

Two things limit how far this can be pushed here.  The model tracks only
photons energetic enough to ionize the gas, because that is what the
photoionization data describes; photons below the ionization threshold but
above the cathode work function would also eject electrons, and they are
simply absent.  And :math:`\vec{\gamma}_\Psi`, like its ion counterpart, is
a surface property carrying the same order-of-magnitude uncertainty.  The
photoemission channel can therefore be scaled with ``xi_emit``, the
photoionization channel with ``xi_photo``, and both switched off entirely
with ``cone_angle: 0`` (:ref:`Chap:Configuration`) — which is how a
calculation can be made to show what the photons are actually contributing.

.. admonition:: Code

   ``get_gamma_Psi(EN, p, T)`` returns :math:`\vec{\gamma}_\Psi`, shape
   :math:`(N_\gamma,)`.  :class:`incept1d.mechanism.Mechanism` multiplies it by
   ``xi_emit`` and the output of ``get_B`` by ``xi_photo``.
