.. _Chap:NewMechanisms:

Creating a new mechanism
========================

.. contents:: On this page
   :local:
   :depth: 1

A *mechanism file* is a plain Python module that describes one gas: its
tracked species, reactions, transport coefficients, photoionization data,
and cathode yields.  It is not imported as a package — :func:`incept1d.mechanism.load_mechanism`
executes it with ``importlib`` after injecting configuration variables into
its namespace — so it must be self-contained and must locate its own data
files relative to ``__file__``.

This page is the interface: what a module must expose, and what each function
must return.  It deliberately does not describe any particular gas.  Two
mechanisms are documented as worked examples — a minimal one in
:ref:`Chap:PaschenMechanism` and a complete one in :ref:`Chap:AirScheme` —
and either is a reasonable starting point to copy.

Required interface
------------------

:func:`incept1d.mechanism.load_mechanism` checks for the following attributes and
refuses to load a module that lacks any of them:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Attribute
     - Contract
   * - ``SPECIES``
     - Ordered list of tracked species names, e.g.
       ``["e", "M+", "M-"]``.  Names are the tokens used in reaction
       strings.  Order defines the row/column order of every matrix.
   * - ``ELECTRON_INDEX``
     - Index of the electron in ``SPECIES``.
   * - ``get_R(EN, p, T, multipliers=None)``
     - Reaction-rate matrix :math:`\bm{R}`, shape :math:`(N_s, N_s)`, in
       s\ :sup:`-1`.  :math:`R_{ij}` is the rate at which one particle of
       species :math:`j` produces (:math:`>0`) or destroys (:math:`<0`)
       species :math:`i`.  Must honour ``multipliers`` (see
       :ref:`Chap:Reactions`).
   * - ``get_V(EN, p, T)``
     - Diagonal drift-velocity matrix :math:`\bm{V}`, shape
       :math:`(N_s, N_s)`, in m/s, with :math:`V_{ii} =
       -\mathrm{sgn}(Z_i)\mu_i|E|` (electrons and anions positive, cations
       negative; see :ref:`Chap:Transport`).  Must be non-singular for all
       :math:`E/N > 0`.
   * - ``get_Pi_e()``, ``get_Pi_plus()``, ``get_Pi_minus()``
     - Row-selection matrices :math:`\bm{\Pi}_\mathrm{e}` (shape
       :math:`(1, N_s)`), :math:`\bm{\Pi}_+` (:math:`(N_+, N_s)`), :math:`\bm{\Pi}_-`
       (:math:`(N_-, N_s)`).  Every species must be selected by exactly one
       of them.
   * - ``get_gamma_plus(EN, p, T)``
     - Ion-induced SEE yields :math:`\vec{\gamma}_+`, shape
       :math:`(N_+,)`, in the row order of ``get_Pi_plus()``.  ``EN`` is the
       cathode field.
   * - ``get_gamma_plus_with(EN, p, T, gamma0=None, gamma1=None, eref=None, beta=None)``
     - Same, with optional per-call overrides (``None`` = module default).
       This is the function the solver actually calls.
   * - ``get_B(EN, p, T)``
     - Photon absorption coupling :math:`\bm{B}`, shape
       :math:`(N_s, N_\gamma)`, :math:`B_{ij} = \beta_j[i]\xi_j\kappa_j` in
       m\ :sup:`-1`.
   * - ``get_C(EN, p, T)``
     - Photon source coupling, shape :math:`(N_\gamma, N_s)`, acting on
       *densities*: :math:`C_{j,\mathrm{e}} = \frac{\Delta\Omega}{4\pi}\rho\,g_j`
       in s\ :sup:`-1`, zero elsewhere.  The solver multiplies by
       :math:`\bm{V}^{-1}`.
   * - ``get_kappa(p, T)``
     - Absorption coefficients :math:`\vec{\kappa}`, shape
       :math:`(N_\gamma,)`, in m\ :sup:`-1`.
   * - ``get_gamma_Psi(EN, p, T)``
     - Photon-induced SEE yields :math:`\vec{\gamma}_\Psi`, shape
       :math:`(N_\gamma,)`.

A mechanism without photoionization returns arrays with :math:`N_\gamma =
0` from ``get_B`` (shape ``(N_s, 0)``), ``get_C`` (``(0, N_s)``),
``get_kappa`` and ``get_gamma_Psi`` (``(0,)``).

Optional interface
------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Attribute
     - Used by
   * - ``alpha(EN, p, T)``, ``eta(EN, p, T)``
     - Townsend coefficients in m\ :sup:`-1`.  Required for
       ``--plot-ionization-integral``, ``--streamer-criterion`` and
       ``incept1d ionization``.
   * - ``init_photoionization(ngroups, cone_angle_deg)``
     - Called by ``Config.post_exec_init`` to (re)build the photon groups.
   * - ``REACTIONS``
     - The declarative reaction list; read by
       ``incept1d chombo`` to export raw rate coefficients.
   * - ``ElectronMeanEnergy(EN)``, ``ElectronMobility(EN)``,
       ``ElectronDiffusion(EN)``, ``<Ion>Mobility(EN)`` …
     - Transport helpers used by ``incept1d chombo``.

Configuration hooks
-------------------

Parameters that must be known *before* the module body runs — which data file
to read, the cathode constants — are read from the module namespace with a
default, so that ``Config.pre_exec_vars()`` can inject them:

.. code-block:: python

   # Read at module level, so an injected value wins over the default.
   CROSS_SECTIONS = globals().get(
       "CROSS_SECTIONS", os.path.join(_HERE, "default_swarm_data.txt")
   )
   _GAMMA0 = globals().get("_GAMMA0", 1.0e-2)

Anything a mechanism computes at import — reading a table, fitting a photon
decomposition — must therefore come *after* these reads.

Parameters that act at run time — reaction multipliers, ``xi_photo``,
``xi_emit``, polarity overrides — are applied by
:class:`incept1d.mechanism.Mechanism` around the interface functions; the mechanism
file does not need to know about them beyond honouring the ``multipliers``
argument of ``get_R``.

A new mechanism *family* (a new directory) needs its own ``config.py``
implementing ``pre_exec_vars`` / ``post_exec_init`` / ``mechanism_params`` /
``label``; :ref:`Chap:ConfigPy` specifies that protocol and gives a minimal
implementation.  The solvers never import it directly.

Step by step
------------

1. **Species.**  Decide the tracked species and their order; put the
   electron first by convention.  Classify each as electron, cation or
   anion — this fixes ``get_Pi_*``.
2. **Swarm data.**  Obtain electron mobility, diffusion, mean energy and the
   impact-ionization and attachment rate coefficients against :math:`E/N`
   for the mixture, usually from a Boltzmann solver, and write a loader for
   whatever format they arrive in.  Interpolate in :math:`\log E/N`.
3. **Ion transport.**  Reduced mobilities :math:`\mu_iN` for every ion,
   constant or tabulated; assemble ``get_V``.
4. **Reactions.**  Write ``REACTIONS`` as ``(string, rate)`` pairs, fold
   neutral densities into the rate, and build ``get_R`` with
   :func:`incept1d.reactions.compile_reactions` /
   :func:`incept1d.reactions.build_R_from_compiled`.
   Remember: exactly one tracked species on the left of every reaction.
5. **Photoionization.**  Either implement :math:`\bm{B}`, :math:`\bm{C}` and
   :math:`\vec{\kappa}` for the gas, or return zero-width arrays to declare
   that the mechanism has none (:ref:`Chap:Photoionization`).
6. **Cathode.**  Implement ``get_gamma_plus_with`` and ``get_gamma_Psi``.
7. **Optional helpers.**  ``alpha``/``eta`` if you want ionization
   integrals and the streamer criterion.
8. **Test.**  Run the standalone plot (``python3 path/to/mechanism.py``),
   then ``incept1d eigenvalues``, then ``incept1d pdiv`` with a configuration that
   reduces the chemistry to the textbook limit and check it against
   :eq:`eq_standard_paschen`.

Common pitfalls
---------------

* **Units.**  ``EN`` is in Td everywhere; ``p`` in bar; ``T`` in K; all
  returns in SI.  The neutral density is :math:`N = p\cdot10^5/(k_BT)`.
* **Sign of** :math:`\bm{V}`.  A wrong sign on an ion velocity turns a
  boundary condition into a nonsense constraint and usually produces no
  roots at all.
* **Zero velocities.**  :math:`\bm{V}^{-1}` is formed by inverting the
  diagonal, so every species needs a non-zero mobility.
* **Import-time work.**  Loading the tables and fitting the photon groups
  happens at import; keep it fast, because every configuration re-executes
  the module.
* **Relative paths.**  Use ``os.path.join(_HERE, ...)`` for data files;
  the working directory is the repository root, not the mechanism
  directory.
