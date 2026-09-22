.. _Chap:AirScheme:

A minimal scheme for dry air
============================

.. contents:: On this page
   :local:
   :depth: 1

The mechanism shipped as ``mechanisms/air/air_pancheshnyi.py`` is the reference
dry-air scheme used in the examples.  It tracks electrons and five ion
species,

.. math::

   \vec{n} = \left(n_\mathrm{e},\; n_{\mathrm{N}_2^+},\; n_{\mathrm{O}_2^+},\;
                   n_{\mathrm{O}^-},\; n_{\mathrm{O}_2^-},\; n_{\mathrm{O}_3^-}\right)^\intercal,

in 79 % N\ :sub:`2` / 21 % O\ :sub:`2` at :math:`T = 300` K, with the
reactions listed in :numref:`tab_reactions`.

Reactions
---------

.. _tab_reactions:

.. list-table:: Plasma-chemistry scheme for dry air.  :math:`\mathrm{M}` in
   reaction 8 is either O\ :sub:`2` or N\ :sub:`2`; :math:`T_\mathrm{e}` is
   the electron temperature; :math:`E/N` in Td.
   :header-rows: 1
   :widths: 4 34 50 12

   * - #
     - Reaction
     - Rate
     - Ref.
   * - **Impact ionization**
     -
     -
     -
   * - 1
     - :math:`\mathrm{e} + \mathrm{N}_2 \xrightarrow{k_1} 2\mathrm{e} + \mathrm{N}_2^+`
     - BOLSIG+
     - [Hagelaar2005]_, [Lisbon]_
   * - 2
     - :math:`\mathrm{e} + \mathrm{O}_2 \xrightarrow{k_2} 2\mathrm{e} + \mathrm{O}_2^+`
     - BOLSIG+
     - [Hagelaar2005]_, [Lisbon]_
   * - **Attachment and detachment**
     -
     -
     -
   * - 3
     - :math:`\mathrm{e} + \mathrm{O}_2 \xrightarrow{k_3} \mathrm{O} + \mathrm{O}^-`
     - BOLSIG+
     - [Hagelaar2005]_, [Lisbon]_
   * - 4
     - :math:`\mathrm{e} + \mathrm{O}_2 + \mathrm{O}_2 \xrightarrow{k_4} \mathrm{O}_2^- + \mathrm{O}_2`
     - :math:`1.4\times10^{-41}\,\frac{300}{T_\mathrm{e}}\exp\!\left(-\frac{600}{T}\right)\exp\!\left(\frac{700(T_\mathrm{e}-T)}{T_\mathrm{e}T}\right)` m\ :sup:`6`/s
     - [Kossyi1992]_
   * - 5
     - :math:`\mathrm{O}_2^- + \mathrm{O}_2 \xrightarrow{k_5} \mathrm{e} + 2\mathrm{O}_2`
     - :math:`1.24\times10^{-17}\exp\!\left[-\left(\frac{259}{14.6 + E/N}\right)^2\right]` m\ :sup:`3`/s
     - [Goodson1974]_, [Pancheshnyi2013]_
   * - 6
     - :math:`\mathrm{O}^- + \mathrm{N}_2 \xrightarrow{k_6} \mathrm{e} + \mathrm{N}_2\mathrm{O}`
     - :math:`3.98\times10^{-17}\left(\frac{E/N}{53}\right)^{-2.72}\exp\!\left[-\left(\frac{177}{E/N}\right)^2\right]` m\ :sup:`3`/s
     - [Shuman2023]_, [MalagonRomero2024]_
   * - **Ion conversion**
     -
     -
     -
   * - 7
     - :math:`\mathrm{O}^- + \mathrm{O}_2 \xrightarrow{k_7} \mathrm{O} + \mathrm{O}_2^-`
     - :math:`6.96\times10^{-17}\exp\!\left[-\left(\frac{198}{5.6 + E/N}\right)^2\right]` m\ :sup:`3`/s
     - [Pancheshnyi2013]_
   * - 8
     - :math:`\mathrm{O}^- + \mathrm{O}_2 + \mathrm{M} \xrightarrow{k_8} \mathrm{O}_3^- + \mathrm{M}`
     - :math:`1.1\times10^{-42}\exp\!\left[-\left(\frac{E/N}{65}\right)^2\right]` m\ :sup:`6`/s
     - [Pancheshnyi2013]_
   * - **Photoionization**
     -
     -
     -
   * - 9
     - :math:`\mathrm{e} + \mathrm{N}_2/\mathrm{O}_2 \xrightarrow{k_9} \hbar\omega`
     - :eq:`eq_photogeneration`
     - [Zheleznyak1982]_, [Bourdon2007]_
   * - 10
     - :math:`\hbar\omega + \mathrm{O}_2 \xrightarrow{\gamma_{10}} \mathrm{e} + \mathrm{O}_2^+`
     - :math:`\gamma_{10} = 0.1`
     - [Zheleznyak1982]_
   * - **Secondary electron emission**
     -
     -
     -
   * - 11
     - :math:`\mathrm{N}_2^+ \xrightarrow{\gamma_{11}} \mathrm{e} + \mathrm{N}_2`
     - :math:`\gamma_{11} = 10^{-3}`
     - [Raizer1991]_
   * - 12
     - :math:`\mathrm{O}_2^+ \xrightarrow{\gamma_{12}} \mathrm{e} + \mathrm{O}_2`
     - :math:`\gamma_{12} = 10^{-3}`
     - [Raizer1991]_
   * - 13
     - :math:`\hbar\omega \xrightarrow{\gamma_{13}} \mathrm{e}`
     - :math:`\gamma_{13} = 0.1`
     - [Raizer1991]_

Reactions 1–8 define the electron-ion dynamics and are the ones assembled
into :math:`\bm{R}`; reactions 9–10 enter through the photon blocks
:math:`\bm{B}`, :math:`\bm{C}`, :math:`\bm{\kappa}`
(:ref:`Chap:Photoionization`); reactions 11–13 are the cathode boundary
condition (:ref:`Chap:SecondaryEmission`).  In the mechanism file the
reactions 1–8 appear verbatim as strings, e.g. ``"O- + N2 -> e + N2O"``,
and are parsed into :math:`\bm{R}` by :mod:`incept1d.reactions`
(:ref:`Chap:ModifyingReactions`).

For the ion-conversion and detachment reactions the data of
[Pancheshnyi2013]_ are used.  Several of these rates were re-analysed by
[Hosl2017]_ with pulsed Townsend experiments; reference data from
different sources do not agree quantitatively, although they follow the
same trends.  The associative-detachment rate (reaction 6) is taken from
[MalagonRomero2024]_ rather than the older compilation, since the original
experimental determination was recently found to be inaccurate
[Shuman2023]_.

Electron and ion transport
--------------------------

Electron transport coefficients and the rates :math:`k_1, k_2, k_3` are
computed with BOLSIG+ [Hagelaar2005]_ from the IST-Lisbon cross sections
[Lisbon]_ by default; the Phelps, Biagi, Trinity, and Morgan sets are
included for sensitivity studies and selected with the ``cross_sections``
key of the JSON configuration.  The Townsend coefficients are

.. math::

   \frac{\alpha}{N} = \frac{\nu_{\mathrm{N}_2}k_1 + \nu_{\mathrm{O}_2}k_2}{\mu_\mathrm{e}E},
   \qquad
   \frac{\eta}{N} = \frac{\nu_{\mathrm{O}_2}\left(k_3 + \nu_{\mathrm{O}_2}k_4N\right)}{\mu_\mathrm{e}E}.

Positive ions use a constant reduced mobility
:math:`\mu_\mathrm{ion}N = 5\times10^{21}\,\mathrm{m^{-1}\,V^{-1}\,s^{-1}}`
[Bohringer1987]_;
:math:`\mathrm{O}^-` uses :math:`1.2\times10^{22}`; the
:math:`\mathrm{O}_2^-` and :math:`\mathrm{O}_3^-` mobilities are tabulated
against :math:`E/N` from LXCat (``mechanisms/air/o2m_mobility.txt``,
``mechanisms/air/o3m_mobility.txt``).

Three-body attachment and the Bloch-Bradbury mechanism
------------------------------------------------------

Reaction 4 is written as a three-body process, but
:math:`\mathrm{O}_2^-` formation is really a sequence of two-body steps —
the Bloch-Bradbury mechanism:

.. math::

   \begin{aligned}
   \mathrm{e} + \mathrm{O}_2 &\xrightarrow{k_c} \mathrm{O}_2^{-*}, \\
   \mathrm{O}_2^{-*} &\xrightarrow{k_\tau} \mathrm{e} + \mathrm{O}_2, \\
   \mathrm{O}_2^{-*} + \mathrm{O}_2 &\xrightarrow{k_d} \mathrm{e} + 2\mathrm{O}_2, \\
   \mathrm{O}_2^{-*} + \mathrm{O}_2 &\xrightarrow{k_s} \mathrm{O}_2^- + \mathrm{O}_2 ,
   \end{aligned}

where :math:`\mathrm{O}_2^{-*}` is a vibrationally excited ion,
:math:`k_c` the capture rate, :math:`k_\tau` the autodetachment rate,
:math:`k_d` the collisional detachment rate, and :math:`k_s` the collisional
stabilisation rate.  Eliminating the excited state gives an effective
three-body rate

.. math::

   k_a = \frac{k_ck_s}{k_\tau + (k_s + k_d)[\mathrm{O}_2]} .

When autodetachment dominates, :math:`k_\tau \gg (k_s + k_d)[\mathrm{O}_2]`,
the dynamics are genuinely three-body with :math:`k_a \approx k_ck_s/k_\tau`
— this is the [Kossyi1992]_ rate used in reaction 4.  When collisional stabilisation dominates the
process becomes an effective two-body attachment with rate
:math:`k_ck_s/(k_s + k_d)`.  With :math:`(k_s + k_d) \sim 10^{-15}`
m\ :sup:`3`/s and :math:`k_\tau \sim 10^{10}` s\ :sup:`-1`, the transition
occurs around 1 bar, at the bottom of the pressure range where the scheme
has been validated, so the three-body form likely *overestimates*
:math:`\mathrm{O}_2^-` formation at high pressure.  A similar caveat applies
to :math:`\mathrm{O}_3^-` formation (reaction 8).

.. admonition:: Code

   ``mechanisms/air/air_2body.py`` is an alternative mechanism that keeps
   :math:`\mathrm{O}_2^{-*}` and :math:`\mathrm{O}_3^{-*}` as explicit
   tracked species (eight species in total) and represents the
   Bloch-Bradbury sequence with two-body reactions.  It shares the same
   interface and configuration protocol and can be substituted for
   ``air_pancheshnyi.py`` in every command.

Stability of :math:`\mathrm{O}_3^-`
-----------------------------------

Reaction 8 is the only reaction in the scheme that acts as a *permanent*
electron sink.  Collisional detachment :math:`\mathrm{O}_3^- \to
\mathrm{O}_3 + \mathrm{e}` has an energy barrier of about 2.1 eV, and
dissociation :math:`\mathrm{O}_3^- \to \mathrm{O}^- + \mathrm{O}_2` about
1.6–1.8 eV, both far above typical negative-ion mean energies of
0.1–0.2 eV at 100 Td.  :math:`\mathrm{O}_3^-` may therefore be stable only
in a moderate sense and begin to contribute to the attachment-detachment
cycle at very long gaps; no experimental rates are available.  The
configuration set ``mechanisms/air/ionsensitivity.json`` switches reactions 7 and 8
off individually to quantify their role.

Pressure scaling
----------------

For a given :math:`E/N` the negative-ion drift velocity is essentially
pressure-independent, so the transit time :math:`d/v_-` is too.  The
collisional detachment rates (reactions 5, 6) scale with :math:`N \propto
p`, so the mean number of detachment events during transit is proportional
to :math:`pd`.  As :math:`pd` grows — through higher pressure (shorter
detachment time) or a longer gap (longer transit) — negative ions are
increasingly likely to return their electrons to the swarm, and those
electrons feed further impact ionization.  This is the physical origin of
the departure from the classical Paschen law at large :math:`pd`.
