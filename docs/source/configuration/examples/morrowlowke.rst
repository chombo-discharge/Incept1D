.. _Chap:MorrowLowkeScheme:

Example: The Morrow–Lowke model for air
=======================================

.. contents:: On this page
   :local:
   :depth: 1

``mechanisms/air/morrowlowke/air_morrowlowke.py`` implements the air model of
[Morrow1997]_, which is widely used in streamer simulations.  It sits between
the two other examples in this chapter.  Like :ref:`Chap:PaschenMechanism` it
needs no data files, because every transport coefficient is an analytic fit
in :math:`E/N`.  Like :ref:`Chap:AirScheme` it has attachment and multigroup
photoionization.  It tracks three species,

.. math::

   \vec{n} = \left(n_\mathrm{e},\; n_+,\; n_-\right)^\intercal,

one electron density and a single lumped positive and negative ion.  The
model has **no detachment**, so a negative ion is a permanent electron loss.

Reactions
---------

.. _tab_reactions_morrowlowke:

.. list-table:: Reaction scheme of the Morrow–Lowke model.  :math:`\mathrm{M}`
   is a neutral molecule; every rate is a Townsend coefficient times the
   electron drift speed :math:`|W_\mathrm{e}|`.
   :header-rows: 1
   :widths: 4 34 38 24

   * - #
     - Reaction
     - Rate (s\ :sup:`-1`)
     - Ref.
   * - 1
     - :math:`\mathrm{e} + \mathrm{M} \to 2\mathrm{e} + \mathrm{M}^+`
     - :math:`\alpha\,|W_\mathrm{e}|`
     - [Morrow1997]_
   * - 2
     - :math:`\mathrm{e} + \mathrm{M} \to \mathrm{M}^-`
     - :math:`\eta_2\,|W_\mathrm{e}|`
     - [Morrow1997]_
   * - 3
     - :math:`\mathrm{e} + 2\mathrm{M} \to \mathrm{M}^- + \mathrm{M}`
     - :math:`\eta_3\,|W_\mathrm{e}|`
     - [Morrow1997]_
   * - 4
     - :math:`\mathrm{e} + \mathrm{N}_2/\mathrm{O}_2 \to \hbar\omega`
     - :eq:`eq_photogeneration`, with :math:`\alpha|W_\mathrm{e}|` as the
       ionization frequency
     - [Zheleznyak1982]_, [Bourdon2007]_
   * - 5
     - :math:`\hbar\omega + \mathrm{O}_2 \to \mathrm{e} + \mathrm{M}^+`
     - :math:`\gamma_5 = 0.1`
     - [Zheleznyak1982]_
   * - 6
     - :math:`\mathrm{M}^+ \to \mathrm{e}` at the cathode
     - :math:`\gamma_6 = 10^{-3}`
     - [Raizer1991]_
   * - 7
     - :math:`\hbar\omega \to \mathrm{e}` at the cathode
     - :math:`\gamma_7 = 0.1`
     - [Raizer1991]_

Reactions 1–3 are assembled into :math:`\bm{R}` from the strings
``"e + M -> 2e + M+"``, ``"e + M -> M-"`` and ``"e + 2M -> M- + M"``, which
are also the keys for ``reaction_multipliers``.  The published model also
has electron-ion and ion-ion recombination and electron diffusion.
Recombination is quadratic in the charged densities and diffusion is not
part of the drift-reaction model (:ref:`Chap:TheoryOverview`), so neither
enters the inception criterion.

Photoionization (reactions 4–5) and secondary emission (reactions 6–7) are
the same as in :ref:`Chap:AirScheme`, including the two-stream fit to the
Zheleznyak absorption curve, the configuration keys and the defaults.  The
only difference is the photon source: it is proportional to the
Morrow–Lowke ionization frequency :math:`\alpha|W_\mathrm{e}|` in place of
the BOLSIG+ rates :math:`k_1, k_2`.

Transport data
--------------

The fits are published with :math:`E/N` in V cm\ :sup:`2`, :math:`N` in
cm\ :sup:`-3` and velocities in cm/s; they are quoted here in those units.
The mechanism converts them to SI and to :math:`E/N` in Td
(1 Td = 10\ :sup:`-17` V cm\ :sup:`2`).

.. list-table:: Morrow–Lowke fits, :math:`x = E/N` in V cm\ :sup:`2`.
   :header-rows: 1
   :widths: 26 44 30

   * - Quantity
     - Fit
     - Range
   * - :math:`|W_\mathrm{e}|` (cm/s)
     - :math:`7.4\times10^{21}x + 7.1\times10^{6}`
     - :math:`x > 2\times10^{-15}`
   * -
     - :math:`1.03\times10^{22}x + 1.3\times10^{6}`
     - :math:`10^{-16} < x \le 2\times10^{-15}`
   * -
     - :math:`7.2973\times10^{21}x + 1.63\times10^{6}`
     - :math:`2.6\times10^{-17} < x \le 10^{-16}`
   * -
     - :math:`6.87\times10^{22}x + 3.38\times10^{4}`
     - :math:`x \le 2.6\times10^{-17}`
   * - :math:`\alpha/N` (cm\ :sup:`2`)
     - :math:`2\times10^{-16}\exp(-7.248\times10^{-15}/x)`
     - :math:`x > 1.5\times10^{-15}`
   * -
     - :math:`6.619\times10^{-17}\exp(-5.593\times10^{-15}/x)`
     - :math:`x \le 1.5\times10^{-15}`
   * - :math:`\eta_2/N` (cm\ :sup:`2`)
     - :math:`8.889\times10^{-5}x + 2.567\times10^{-19}`
     - :math:`x > 1.05\times10^{-15}`
   * -
     - :math:`6.089\times10^{-4}x - 2.893\times10^{-19}`
     - :math:`x \le 1.05\times10^{-15}`
   * - :math:`\eta_3/N^2` (cm\ :sup:`5`)
     - :math:`4.7778\times10^{-59}x^{-1.2749}`
     - All :math:`x`
   * - :math:`\mu_+` (cm\ :sup:`2`/V s)
     - :math:`2.34`
     - All :math:`x`
   * - :math:`\mu_-` (cm\ :sup:`2`/V s)
     - :math:`2.7`
     - :math:`x > 5\times10^{-16}`
   * -
     - :math:`1.86`
     - :math:`x \le 5\times10^{-16}`

Three details of the implementation:

* The ion mobilities are published for atmospheric density.  They are
  taken to hold at :math:`N_0 = 2.5\times10^{19}` cm\ :sup:`-3` and scaled
  as :math:`\mu N = \mathrm{const}`, so every drift velocity depends on
  :math:`E/N` alone.
* The low-field branch of :math:`\eta_2/N` turns negative below
  :math:`x \approx 4.75\times10^{-16}` (47.5 Td).  It is clipped at zero.
* The fits are used as published, including two steps: :math:`|W_\mathrm{e}|`
  drops by about 1.3 % across :math:`x = 10^{-16}` (10 Td), and
  :math:`\mu_-` jumps from 1.86 to 2.7 at 50 Td.  The other breakpoints
  are continuous to better than 0.3 %.

Because :math:`\eta_3 \propto N^2`, the critical field where
:math:`\alpha = \eta` depends on pressure: about 108 Td at 1 bar and 293 K,
rising with pressure.

Checks
------

The model has no detachment, so with photon feedback switched off
(``"cone_angle": 0``) the reduced model :eq:`eq_standard_paschen` applies
exactly, and the test suite checks the full solve against it.  At higher
:math:`pd` the uniform-field inception curve levels off just above the
critical field, since no detachment feedback can pull it below.  In
:ref:`Chap:AirScheme`, detachment does.

The same commands as for the other air schemes apply:

* ``python3 mechanisms/air/morrowlowke/air_morrowlowke.py`` — Plots :math:`\alpha`,
  :math:`\eta_2`, :math:`\eta_3` and :math:`|W_\mathrm{e}|` against
  :math:`E/N`.
* ``incept1d eigenvalues mechanisms/air/morrowlowke/air_morrowlowke.py`` — The leading
  eigenvalue of :math:`\bm{R}\bm{V}^{-1}` against :math:`E/N`.
* ``incept1d pdiv mechanisms/air/morrowlowke/air_morrowlowke.py`` — The inception curve.

``mechanisms/air/morrowlowke/config.py`` is the shared air configuration
class without the ``cross_sections`` key: the transport is analytic, so
there are no swarm tables to select, and a configuration that names one was
written for another mechanism and is rejected.

``incept1d chombo`` does not accept this mechanism: its tables include the
mean electron energy, and the Morrow–Lowke model does not provide one.
