.. _Chap:MechanismFiles:

Mechanism files
===============

.. contents:: On this page
   :local:
   :depth: 1

A mechanism file is a JSON file that describes a gas: its tracked species,
their charges and mobilities, the reactions between them with their rate
coefficients, the Townsend coefficients, the photoionization model and the
electrode yields.  ``Incept1D`` compiles it into the reaction matrix
:math:`\bm{R}`, the drift velocities :math:`\bm{V}`, the photon couplings
:math:`\bm{B}`, :math:`\bm{C}` and :math:`\bm{\kappa}` and the boundary
conditions of :ref:`Chap:SecondaryEmission` — everything a Python mechanism
(:ref:`Chap:NewMechanisms`) builds by hand follows from the file.

Every coefficient is a formula.  The file is checked when it is loaded,
and anything that is not a mechanism — an unknown key, a misspelt name in a
formula, a reaction that does not conserve charge — is an error then,
before a single field is solved.

A first example
---------------

An electronegative gas of two components :math:`A` and :math:`B`, with
ionization of :math:`A`, attachment to :math:`B` and detachment:

.. code-block:: json

   {
       "format": "incept1d-mechanism/1",
       "description": "Two-component electronegative gas, for illustration",
       "parameters": {
           "x_b": 0.2,
           "k_ion0": 3e-15,
           "k_att": 2e-17,
           "gamma0": 1e-3
       },
       "gas": {"A": "1 - x_b", "B": "x_b"},
       "definitions": {
           "k_ion": "k_ion0*exp(-600/EN)",
           "v_e": "4e23*EN*1e-21"
       },
       "alpha": "k_ion*x_A*N/v_e",
       "eta": "k_att*x_B*N/v_e",
       "species": [
           {"name": "e",  "charge": -1, "mobility": "4e23"},
           {"name": "A+", "charge":  1, "mobility": "2e21"},
           {"name": "B-", "charge": -1, "mobility": "2e21"}
       ],
       "reactions": [
           {"reaction": "e + A -> 2e + A+", "rate": "k_ion"},
           {"reaction": "e + B -> B-",      "rate": "k_att"},
           {"reaction": "B- + M -> e + B + M", "rate": "1e-18*exp(-(100/EN)**2)"}
       ],
       "emission": {"ions": {"A+": "gamma0"}}
   }

Run it as any mechanism, ``incept1d pdiv gas.json``, optionally with
configuration files that set its parameters (:ref:`Chap:Configuration`).

Structure
---------

.. list-table::
   :header-rows: 1
   :widths: 20 12 68

   * - Key
     - Required
     - Content
   * - ``format``
     - yes
     - ``"incept1d-mechanism/1"``.  A JSON file without it is not a
       mechanism, which is how a mechanism file is told from a configuration
       file.
   * - ``description``
     - no
     - Free text.
   * - ``parameters``
     - no
     - Named numbers and file names: the defaults of everything a
       configuration may change.
   * - ``gas``
     - no
     - Mole fraction of each neutral component.
   * - ``tables``
     - no
     - Swarm data read from files, each column a variable.
   * - ``untracked``
     - no
     - Products of reactions that are not followed.
   * - ``definitions``
     - no
     - Named formulas, in order.
   * - ``alpha``, ``eta``
     - no
     - The Townsend ionization and attachment coefficients [1/m].
   * - ``species``
     - yes
     - The tracked species, with charge and mobility.
   * - ``reactions``
     - no
     - Reaction strings with rate coefficients.
   * - ``photoionization``
     - no
     - The photon groups and how they couple to the species.
   * - ``emission``
     - no
     - The ion-induced secondary emission yields of the cathode.

Keys beginning with ``_`` are comments and are ignored at every level, so
``"_k_att": "dissociative attachment, m^3/s"`` next to ``"k_att"`` documents
it.  Any other key that is not in the format is an error.

Formulas
--------

A coefficient is a formula in a string, ``"1.24e-17*exp(-(259/(14.6 + EN))**2)"``,
or a plain number.  A formula may use numbers, names, the operators
``+ - * / **``, parentheses, and the functions ``exp``, ``log``, ``log10``,
``sqrt``, ``abs``, ``tanh``, ``cos``, ``sin`` (one argument each), ``min``
and ``max`` (two).  Nothing else is allowed — no comparisons, conditionals,
attributes or other calls — so a formula computes a number and cannot do
anything else; a threshold is written with ``max`` or ``min``.  Arithmetic
follows IEEE floating point: a division by zero gives infinity and
:math:`e^{-\infty}` gives zero, rather than an error, while an undefined
result such as :math:`0/0` is NaN and is warned about.

The names a formula may use are, in order of definition:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Name
     - Value
   * - ``EN``
     - Reduced field :math:`E/N` [Td]
   * - ``E``
     - Field [V/m]
   * - ``N``
     - Neutral number density [m\ :sup:`-3`], from the pressure and
       temperature of the run
   * - ``p``, ``T``
     - Pressure [bar] and gas temperature [K]
   * - ``kB``, ``qe``, ``pi``
     - Boltzmann constant [J/K], elementary charge [C], :math:`\pi`
   * - parameters
     - Every numeric parameter
   * - ``x_<X>``
     - The mole fraction of gas component ``X``
   * - table columns
     - Every column declared in ``tables``, at the current :math:`E/N`
   * - definitions
     - Every definition before the one being defined
   * - ``alpha``, ``eta``
     - When the file gives them: in species, reactions, photoionization and
       emission, not in definitions

Parameters
----------

.. code-block:: json

   "parameters": {
       "gamma0": 1e-3,
       "cross_sections": "data/swarm.txt"
   }

A parameter is a number or a file name.  Parameters are what a configuration
can change, and the only things it can change: a configuration that names
anything else is refused.  So a coefficient that a study will vary belongs
here, with its default, rather than as a literal inside a formula.

A file name is relative to the mechanism file.  When a configuration sets
it, a relative name is looked for beside the configuration file first, and
then beside the mechanism file.

Names start with a letter and contain letters, digits and ``_``; they may
not be one of the built-in names or functions.

Gas
---

.. code-block:: json

   "gas": {"A": "1 - x_b", "B": "x_b"}

The mole fraction of each neutral component, as a formula of the
parameters.  The fractions must lie in :math:`[0, 1]` and sum to 1 to within
:math:`10^{-9}`.  Each
defines ``x_<component>`` for the formulas, and the component names are
what reactions use for neutral partners.

Tables
------

.. code-block:: json

   "tables": {
       "swarm": {
           "format": "bolsig+",
           "file": "cross_sections",
           "columns": {"mean_energy": "A1", "mu_e_N": "A2"},
           "rates": {"k_ion": "A Ionization", "k_att": "B Attachment"}
       },
       "ions": {"format": "lxcat", "file": "data/ion_mobility.txt", "columns": {"K0": 1}}
   }

Swarm data tabulated against :math:`E/N`.  ``file`` is a file name or the
name of a file parameter.  Each entry of ``columns`` and ``rates`` defines a
variable: the column at the :math:`E/N` of the evaluation, interpolated
linearly and held at its end values outside the table.

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Format
     - Content
   * - ``bolsig+``
     - BOLSIG+ output.  ``columns`` names transport coefficients by their
       label in the output (``A1`` mean energy [eV], ``A2`` mobility
       :math:`\times N` [1/(m V s)], ...).  ``rates`` names rate coefficients
       [m\ :sup:`3`/s] by target and process type as the output's legend
       lists them, ``"A Ionization"``; every process of that kind is summed.
   * - ``lxcat``
     - A two-column LXCat table between its first pair of ``----`` lines,
       such as tabulated ion mobilities; column ``1`` is the value.

A BOLSIG+ table is computed for one gas mixture and temperature.  When the
file records its mixture, ``gas`` must agree with it to :math:`10^{-3}`,
since the neutral densities in the reactions must be those the rates were
computed for.  When it does not, the mixture cannot be checked, and a
warning says so unless the configuration sets the mixture itself.  A run
whose temperature differs from the table's by more than 20 K is warned
about once.

Definitions
-----------

.. code-block:: json

   "definitions": {
       "Te": "2*mean_energy*qe/(3*kB)",
       "k_3b": "1.4e-41*(300/Te)"
   }

Named intermediate formulas, evaluated in order at every
:math:`(E/N, p, T)`: a definition may use those before it.  A definition is
a convenience, not a parameter — a configuration cannot change it.

Species
-------

.. code-block:: json

   "species": [
       {"name": "e",  "charge": -1, "mobility": "mu_e_N"},
       {"name": "A+", "charge":  1, "mobility": "2e21"}
   ]

Each tracked species has a name, a non-zero integer charge and a reduced
mobility :math:`\mu N` [1/(m V s)].  The electron ``e``, with charge
:math:`-1`, is required.  The drift velocity is

.. math::

   V_{ii} = -\operatorname{sgn}(Z_i)\, (\mu N)_i\, E/N,

so electrons and negative ions drift toward the anode and positive ions
toward the cathode, and the row selectors :math:`\bm{\Pi}_e`,
:math:`\bm{\Pi}_+` and :math:`\bm{\Pi}_-` of :ref:`Chap:SecondaryEmission`
follow from the charges.  A neutral species does not drift and cannot be
tracked.

Townsend coefficients
---------------------

``alpha`` and ``eta`` [1/m] are given, not derived from the reactions.
They are what ``incept1d ionization`` and ``--streamer-criterion`` use; the
inception criterion itself uses the reactions.  Deriving them would assume
that every ionization is one reaction ``e + X -> 2e + X+``, which is not true
of a model where ionization runs through an intermediate state.  A
reaction multiplier does not change them.  A mechanism without them is
valid; those two features are then unavailable.

Reactions
---------

.. code-block:: json

   "reactions": [
       {"reaction": "e + A -> 2e + A+",        "rate": "k_ion"},
       {"reaction": "e + 2B -> B- + B",        "rate": "k_3b"},
       {"reaction": "B- + B + M -> B3- + M",   "rate": "1e-42"}
   ]

A reaction string has one tracked species on the left, once; the other
reactants are gas components or ``M``, which stands for the total neutral
density.  Tokens are separated by ``" + "``, so ``A+`` is a name, and an
integer prefix is a stoichiometric coefficient, ``2B``.  A product is a
species, a gas component, ``M``, or a name listed in ``untracked``.  Charge
must be conserved.

The rate is the rate coefficient :math:`k` in SI units, and the densities
of the neutral reactants are multiplied in: for a reaction driven by species
:math:`j` with neutral reactants :math:`X` of multiplicity :math:`\nu_X`,

.. math::

   K = k \prod_X n_X^{\nu_X}, \qquad n_X = x_X N, \quad n_M = N,

so :math:`k` is in m\ :sup:`3`/s for one neutral partner, m\ :sup:`6`/s for
two, and s\ :sup:`-1` for none.  Each reaction adds
:math:`(\nu^\mathrm{prod}_i - \nu^\mathrm{reac}_i) K` to :math:`R_{ij}`.
The reaction string is the reaction's name: a configuration scales a
reaction by it (:ref:`Chap:ModifyingReactions`), so two reactions may not
have the same string.

Photoionization
---------------

.. code-block:: json

   "photoionization": {
       "absorption": {"model": "zheleznyak", "groups": "ngroups", "absorber": "B"},
       "source": {
           "reactions": ["e + A -> 2e + A+"],
           "efficiency": "0.6*(1 - cos(45*pi/180))/2"
       },
       "products": {"e": "0.1", "B+": "0.1"},
       "cathode yield": "0.1",
       "direction": "both"
   }

The photon groups of the two-stream model
(:eq:`eq_two_stream`).  ``absorption`` gives the number of groups
:math:`N_\gamma` (a formula of the parameters), their absorption
coefficients and fractions :math:`g_j` (:math:`\sum_j g_j = 1`) from an
absorption model, and the gas component that absorbs.  The one model is
``zheleznyak``, a fit of :math:`N_\gamma` exponentials to the absorption
function of Zheleznyak et al. (:ref:`Chap:AirScheme` shows the fit), with
coefficients :math:`\hat\kappa_j` per unit absorber pressure [1/(m Pa)],
so that, with :math:`p` in Pa,

.. math::

   \kappa_j = \hat\kappa_j\, x_\mathrm{abs}\, p.

The photons are emitted in proportion to the source reactions:

.. math::

   C_{ji} = \varepsilon \, g_j \sum_{r \to i} K_r,

summed over the source reactions driven by species :math:`i`, with
:math:`\varepsilon` the ``efficiency`` formula.  :math:`K_r` is the rate
before reaction multipliers.  ``products`` gives the species an absorbed
photon makes and how many, so that :math:`B_{ij} = c_i \kappa_j`;
``cathode yield`` is the photoemission yield :math:`\gamma_{\Psi,j}` of
every group; ``direction``, ``"both"`` by default, is the
``PHOTON_DIRECTIONS`` of every group (:eq:`eq_two_stream_directional`).

Emission
--------

.. code-block:: json

   "emission": {"ions": {"A+": "gamma0 + gamma1*exp(-eref/(beta*E))"}}

The secondary-emission yield :math:`\gamma_+` of each positive species at
the cathode.  Every positive species needs one; write ``"0"`` for none, so
that a forgotten yield is not mistaken for a zero one.

What is checked when a file is loaded
-------------------------------------

* **Keys.**  Every key is one the format defines, or a ``_`` comment.
* **Names.**  Every name in a formula is defined before it is used; names
  are not reused.
* **Reactions.**  One tracked reactant; known neutral partners and products;
  charge conserved; no two with the same string.
* **Species.**  Unique names, non-zero integer charges, the electron present.
* **Gas.**  Mole fractions in :math:`[0, 1]` summing to 1, and equal to a
  BOLSIG+ table's recorded mixture.
* **Tables.**  The file exists, and every column and rate it is asked for.
* **Photoionization.**  A known model, a whole number of groups, an absorber
  in the gas, source reactions that exist.
* **Emission.**  A yield for every positive species.
* **Configurations.**  Only parameters, ``label``, ``reaction_multipliers``,
  ``positive`` and ``negative``; numbers for numbers and file names for
  files; multipliers that name reactions of this file.

When the format is not enough
-----------------------------

A mechanism that needs what a formula cannot say — a rate that switches
form at a threshold, a coupling that is not the product of a source and an
absorption, a model that is not a reaction list at all — is written in
Python instead (:ref:`Chap:NewMechanisms`).  The solver sees both the same
way.

API reference
-------------

.. automodule:: incept1d.declarative
   :members: load_json_mechanism, is_mechanism_file, DeclarativeMechanism, MechanismFileError, FORMAT

.. automodule:: incept1d.expressions
   :members: Expression, ExpressionError, FUNCTIONS

.. automodule:: incept1d.tables
   :members: read_table, Table, TableError

.. automodule:: incept1d.zheleznyak
   :members: zheleznyak, fit_twostream, twostream_groups
