.. _Chap:RFDScheme:

Example: Relativistic feedback discharge
========================================

.. contents:: On this page
   :local:
   :depth: 1

``mechanisms/air/rfd/`` implements the relativistic fluid model of
[Faerder2026b]_ for the relativistic feedback discharge (RFD) in the
high-field region of a thundercloud.  Nothing in it is a gas-discharge
species in the usual sense, yet it fits the interface unchanged, and it is
the shipped example of two features no other mechanism uses: photon groups
that radiate in one direction only (:eq:`eq_two_stream_directional`), and
feedback that closes entirely in the gas, with nothing emitted at either end.

.. code-block:: text

   mechanisms/air/rfd/
   ├── air_rfd.py     the mechanism module
   ├── config.py      attenuation lengths and the sea-level reference density
   └── radius.json    the two radii of the high-field region in the paper

The model
---------

Relativistic runaway electrons (RREs) avalanche along the field and emit
bremsstrahlung photons, forward and backward.  Forward photons make new
RREs by Compton scattering and photoabsorption, and positrons by pair
production; the positrons run back against the field and make RREs by
Bhabha scattering, and so do the backward photons.  With the densities
:math:`n_\mathrm{r}` of RREs, :math:`n_\mathrm{pos}` of positrons and
:math:`n_{\gamma\uparrow}`, :math:`n_{\gamma\downarrow}` of forward and
backward photons, all moving at :math:`\tilde{c} = 0.89\,c`, the model is

.. math::
   :label: eq_rfd_model

   \begin{aligned}
   \frac{1}{\tilde{c}}\partial_t n_\mathrm{r} + \partial_z n_\mathrm{r}
     &= \lambda^{-1} n_\mathrm{r} + P_\mathrm{bb}\,n_\mathrm{pos}
        + P_{\mathrm{cp}\uparrow}\,n_{\gamma\uparrow}
        + P_{\mathrm{cp}\downarrow}\,n_{\gamma\downarrow}, \\
   \frac{1}{\tilde{c}}\partial_t n_\mathrm{pos} - \partial_z n_\mathrm{pos}
     &= P_{\mathrm{pp}\uparrow}\,n_{\gamma\uparrow}
        - \lambda_\mathrm{a,pos}^{-1} n_\mathrm{pos}, \\
   \frac{1}{\tilde{c}}\partial_t n_{\gamma\uparrow} + \partial_z n_{\gamma\uparrow}
     &= P_{\mathrm{br}\uparrow}\,n_\mathrm{r}
        - \lambda_{\mathrm{a}\uparrow}^{-1} n_{\gamma\uparrow}, \\
   \frac{1}{\tilde{c}}\partial_t n_{\gamma\downarrow} - \partial_z n_{\gamma\downarrow}
     &= P_{\mathrm{br}\downarrow}\,n_\mathrm{r}
        - \lambda_{\mathrm{a}\downarrow}^{-1} n_{\gamma\downarrow},
   \end{aligned}

on :math:`0 \le z \le L`, with nothing entering:
:math:`n_\mathrm{r}(0) = n_{\gamma\uparrow}(0) = 0` and
:math:`n_\mathrm{pos}(L) = n_{\gamma\downarrow}(L) = 0`.  The discharge is
self-sustained above the field :math:`E_\mathrm{s}(L)` at which a
non-trivial steady state exists.

How it maps onto the interface
------------------------------

The coordinate is :math:`x = z`: RREs enter at :math:`x = 0` and drift to
:math:`+x`, so the RRE is the mechanism's electron, and positrons, which
drift to :math:`-x` and must vanish at :math:`x = d`, are its positive ion.
Fluxes are :math:`\tilde{c}` times densities.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Interface
     - Value
   * - ``SPECIES``
     - ``["rre", "e+"]``, with :math:`\bm{V} = \diag(\tilde{c}, -\tilde{c})`
   * - ``get_R``
     - :math:`\tilde{c}\begin{bmatrix}\lambda^{-1} & P_\mathrm{bb}\\
       0 & -\lambda_\mathrm{a,pos}^{-1}\end{bmatrix}`
   * - ``PHOTON_DIRECTIONS``
     - ``["forward", "backward"]``: one group per photon population
   * - ``get_kappa``
     - :math:`(\lambda_{\mathrm{a}\uparrow}^{-1},\;
       \lambda_{\mathrm{a}\downarrow}^{-1})`
   * - ``get_C``
     - :math:`\tilde{c}\,P_{\mathrm{br}\uparrow}` and
       :math:`\tilde{c}\,P_{\mathrm{br}\downarrow}` per RRE density
   * - ``get_B``
     - :math:`\begin{bmatrix}P_{\mathrm{cp}\uparrow} &
       P_{\mathrm{cp}\downarrow}\\ P_{\mathrm{pp}\uparrow} & 0\end{bmatrix}`
   * - ``get_gamma_plus``, ``get_gamma_Psi``
     - Zero: nothing is emitted at either end

With these, the steady state of :eq:`eq_rfd_model` is the augmented ODE
:eq:`eq_augmented_ode` term for term, and the open ends are the boundary
conditions of :ref:`Chap:SecondaryEmission` with every yield zero.  The
inception criterion is therefore exactly the self-sustained condition of
the model, and ``incept1d pdiv`` computes :math:`E_\mathrm{s}(L)`.  Because
nothing is emitted at :math:`x = 0`, the Riccati criterion sees the threshold
only as a pole of the reflection operator (:ref:`Chap:Numerics:Riccati`).

Field and density
-----------------

The paper gives every coefficient as a function of the field magnitude
:math:`E` at sea-level density, scaled by the relative density
:math:`\tilde{n}`, with every threshold field scaling as :math:`\tilde{n}`.
The mechanism evaluates each one at the sea-level-equivalent field
:math:`E/\tilde{n}` and multiplies by :math:`\tilde{n}`, which is the paper's
expression at :math:`\tilde{n} = 1` and depends only on :math:`E/N`
elsewhere, so the curve depends only on :math:`N d` like every other
mechanism's.  The reference density :math:`N_0` with
:math:`\tilde{n} = N/N_0` sets the conversion between :math:`E/N` and
:math:`E`; it defaults to 1 atm and 293.15 K and is a configuration key.

Coefficients
------------

At :math:`\tilde{n} = 1`, with :math:`E` in V/m and lengths in m, from the
fits to GEANT4 simulations in [Faerder2026b]_:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Coefficient
     - Value
   * - :math:`\lambda^{-1}`
     - :math:`(E - E_\mathrm{d})/U_\mathrm{R}`, :math:`U_\mathrm{R} =
       7300` kV, :math:`E_\mathrm{d} = 276` kV/m above
       :math:`E_\mathrm{t} = 284` kV/m and 320 kV/m below
   * - :math:`P_\mathrm{bb}`
     - :math:`f(E;\,4.325\times10^{-3},\,7.084\times10^{-10},\,305.4,\,427.6,\,1.1)`
   * - :math:`P_{\mathrm{cp}\uparrow}`
     - :math:`f(E;\,2.468\times10^{-3},\,1.028\times10^{-8},\,-623.9,\,1452.9,\,5.5)`
   * - :math:`P_{\mathrm{cp}\downarrow}`
     - :math:`f(E;\,1.019\times10^{-3},\,1.062\times10^{-8},\,-10.7,\,874.9,\,4.1)`
   * - :math:`P_{\mathrm{pp}\uparrow}`
     - :math:`8.580\times10^{-5}`
   * - :math:`P_{\mathrm{br}\uparrow}`, :math:`P_{\mathrm{br}\downarrow}`
     - :math:`1.442\times10^{-2}`, :math:`9.666\times10^{-4}`
   * - :math:`\lambda_\mathrm{a,pos}^{-1}`
     - :math:`127/(E - 259.1\times10^{3})`
   * - :math:`\lambda_{\mathrm{a}\uparrow}`, :math:`\lambda_{\mathrm{a}\downarrow}`
     - 307.64, 92.22 for a wide region; 275.31, 55.29 for a radius of 100 m

with

.. math::

   f(E;\,P,\,s,\,E_\mathrm{d},\,E_\mathrm{s},\,r) =
   (P + sE)\left(1 - \exp\left[-\left(\frac{E - E_\mathrm{d}}
   {E_\mathrm{s}}\right)^{r}\right]\right)

for :math:`E > E_\mathrm{d}` and zero below, :math:`E_\mathrm{d}` and
:math:`E_\mathrm{s}` given in kV/m.  The photon attenuation lengths depend on
the horizontal radius :math:`R` of the high-field region, because photons
escape a narrow region sideways; ``radius.json`` holds both cases as two
configurations.

Below 305.4 kV/m there is no Bhabha production, so positrons feed nothing
back and their attenuation does not matter; the mechanism holds
:math:`\lambda_\mathrm{a,pos}^{-1}` at its 305.4 kV/m value there instead of
letting it diverge at 259.1 kV/m.

Normal modes
------------

The steady state of :eq:`eq_rfd_model` is
:math:`\partial_z\vec{n} = \bm{M}\vec{n}` for
:math:`\vec{n} = (n_\mathrm{r}, n_\mathrm{pos}, n_{\gamma\uparrow},
n_{\gamma\downarrow})`, with

.. math::
   :label: eq_rfd_normal_modes

   \bm{M} = \begin{bmatrix}
     \lambda^{-1} & P_\mathrm{bb} & P_{\mathrm{cp}\uparrow} &
     P_{\mathrm{cp}\downarrow} \\
     0 & \lambda_\mathrm{a,pos}^{-1} & -P_{\mathrm{pp}\uparrow} & 0 \\
     P_{\mathrm{br}\uparrow} & 0 & -\lambda_{\mathrm{a}\uparrow}^{-1} & 0 \\
     -P_{\mathrm{br}\downarrow} & 0 & 0 & \lambda_{\mathrm{a}\downarrow}^{-1}
   \end{bmatrix},

whose solutions are sums of modes :math:`\vec{v}_j e^{\kappa_j z}`,
:math:`\bm{M}\vec{v}_j = \kappa_j\vec{v}_j`.  The four open-end conditions
admit a non-trivial combination where the determinant of the modes'
boundary values vanishes; divided by :math:`\det[\vec{v}_1 \cdots
\vec{v}_4]`, that determinant is real and does not depend on how the modes
are scaled or ordered.  ``normal_mode_threshold`` in the mechanism module
evaluates it, and is what the example compares the solver against
(:ref:`Chap:Examples:RFD`).  The :math:`(3,3)` entry is written from the
model equations; it is :math:`-\lambda_{\mathrm{a}\uparrow}^{-1}`, as the
forward photons are attenuated, and so it is in the authors' own
normal-mode code.

Limits of the mapping
---------------------

* **Growth rates.**  :ref:`Chap:Lambda` retards photons with
  :math:`\lambda/c`; the model moves them at :math:`\tilde{c}`, so a growth
  rate above threshold differs from the model's by the photon term.  The
  threshold itself (:math:`\lambda = 0`) is unaffected.
* **Field-dependent coefficients only.**  The paper fits each coefficient in
  a uniform field.  A non-uniform ``--field`` applies them locally, which is
  the model's own assumption but has not been checked against GEANT4.
* **No Townsend coefficients.**  The mechanism defines no ``alpha`` or
  ``eta``, so ``incept1d ionization`` and ``--streamer-criterion`` do not
  apply.
