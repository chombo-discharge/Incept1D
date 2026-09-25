.. _Chap:Numerics:Propagator:

The propagator
==============

.. contents:: On this page
   :local:
   :depth: 1

The formal solution :eq:`eq_theta_soln` requires the propagator
:math:`\bm{M}(d)` from the cathode to the anode.  This page describes how it
is computed.

Uniform field
-------------

For a uniform field :math:`\bm{\mathcal{A}}` is constant and
:math:`\bm{M}(d) = e^{\bm{\mathcal{A}}d}` is evaluated with
:func:`scipy.linalg.expm` (Padé approximation with scaling and squaring).
The full augmented matrix is small — a handful of species and photon groups,
so :math:`N_s + 2N_\gamma` is of order ten — and a single ``expm`` costs
microseconds.

This is the *exact* propagator, not a quadrature of it, so the grid options
have nothing to refine: ``--dx`` and ``--method`` are accepted but ignored
for a uniform field, and every setting returns bit-for-bit the same
:math:`\det\bm{Q}`.  :func:`incept1d.solver.inception_det` takes this
shortcut before building the segment grid, which makes a uniform-field
:math:`pd` sweep roughly an order of magnitude cheaper than the composed
path below (one ``expm`` and one :math:`\bm{\mathcal{A}}` assembly per
:math:`\det\bm{Q}`, against 15 and 16 for the five-segment default).

Non-uniform field: composed propagators
---------------------------------------

For a non-uniform field the domain :math:`[0, d]` is split into
:math:`N_{\min}` initial segments (``--dx N_min``, default 5).  On each
segment :math:`[x_i, x_{i+1}]` of width :math:`h` a local propagator
:math:`\bm{P}_i` is computed, and the total propagator is the path-ordered
product

.. math::

   \bm{M}(d) = \bm{P}_{N}\cdots\bm{P}_2\bm{P}_1 .

Two local propagators are available (``--method``):

**Midpoint rule** (zeroth-order Magnus, default)

.. math::

   \bm{P}_i = \exp\left(h\,\bm{\mathcal{A}}(x_{i+1/2})\right).

It is second-order accurate in :math:`h` and, because ``expm`` handles
arbitrary step sizes, robust for any :math:`h`.  Implemented in
:func:`incept1d.solver.midpoint_propagator`.

**Second-order Magnus** (``magnus2``) with two-point Gauss-Legendre
quadrature, :math:`\bm{\mathcal{A}}_{1,2} = \bm{\mathcal{A}}(x_{i+1/2} \mp
h/(2\sqrt{3}))`:

.. math::

   \bm{\Omega} = \frac{h}{2}\left(\bm{\mathcal{A}}_1 + \bm{\mathcal{A}}_2\right)
               + \frac{\sqrt{3}\,h^2}{12}\left[\bm{\mathcal{A}}_2, \bm{\mathcal{A}}_1\right],
   \qquad \bm{P}_i = e^{\bm{\Omega}} .

The first term is the quadrature approximation to
:math:`\int\bm{\mathcal{A}}\,dx`; the commutator is the leading Magnus
correction and vanishes for a uniform field, so ``magnus2`` reduces exactly
to the midpoint rule there.  It is fourth-order accurate but the Magnus
series only converges when :math:`\int\|\bm{\mathcal{A}}\|\,dx < \pi`; if
:math:`\|\bm{\Omega}_1\|_F > \pi` the code falls back to the midpoint rule
for that segment.  Implemented in :func:`incept1d.solver.magnus2_propagator`.
``magnus2`` uses the fixed :math:`N_{\min}` grid with no adaptation.

Adaptive step halving
---------------------

With the midpoint propagator each initial segment is refined adaptively
(:func:`incept1d.solver._adaptive_midpoint_segment`): the one-step propagator
:math:`\bm{P}_\mathrm{coarse}` over :math:`[x_i, x_{i+1}]` is compared with
the two-step estimate :math:`\bm{P}_\mathrm{fine} = \bm{P}_\mathrm{r}
\bm{P}_\mathrm{l}` over the two halves, and the segment is accepted when

.. math::

   \frac{\|\bm{P}_\mathrm{fine} - \bm{P}_\mathrm{coarse}\|_F}
        {\|\bm{P}_\mathrm{fine}\|_F} \le \mathrm{tol}

(``--dx`` third argument, default 0.03).  Otherwise each half is refined
recursively, re-using the already computed half-step propagators as the
children's coarse estimates so that each level costs only the new
half-steps.  The recursion depth is limited by
:math:`\lfloor\log_2(N_{\max}/N_{\min})\rfloor` (``--dx`` second argument,
default 200), so the total number of fine steps never exceeds
:math:`N_{\max}`.  Setting :math:`N_{\max} = N_{\min}` disables adaptation and
gives a fixed uniform grid.

The step halving judges a segment by comparing one midpoint step with two
half steps, and cannot see a feature that falls between those sample
points: a thin high-field layer at a small electrode can be missed
altogether, and with it the avalanche.  The :math:`N_{\min}` initial
segments are therefore first split until the field changes by at most
20 % across each (:func:`incept1d.solver._field_following_edges`), which
for a wire of radius :math:`a` in a gap of :math:`100\,a` adds a couple of
dozen short segments at the wire; a uniform or gently varying field keeps
its :math:`N_{\min}` equal segments.  With that start the defaults keep
the inception field within about 0.1–0.2 % of its converged value on
sphere-plane and thin-wire coaxial gaps; ``--dx 5 1000 1e-3`` brings that
to about 0.01 % at three to seven times the cost.

Eigenvalue shifting
-------------------

Over a long gap :math:`e^{\bm{\mathcal{A}}d}` has entries of order
:math:`e^{\lambda_{\max}d}`, which overflows double precision well before
:math:`pd` reaches the range of interest.  :func:`incept1d.solver._expm_shifted`
therefore evaluates :math:`\exp(\bm{\Omega} - \mu\bm{I})` with
:math:`\mu = \max(0, \max_i\mathrm{Re}\,\lambda_i(\bm{\Omega}))`.  The
omitted factor :math:`e^{\mu}` is a positive scalar common to every column
of the propagator, hence to every row of :math:`\bm{Q}_d`; it scales
:math:`\det\bm{Q}` by a positive constant and leaves its **sign** — the
only thing the root finder needs — unchanged.  The factor is intentionally
not restored.

Photon groups
-------------

Every photon group is propagated explicitly, however large its optical
depth :math:`\kappa_j d`.  Replacing a thick group by its steady-state
contribution :math:`2\vec{b}_j\vec{c}_j^\intercal/\kappa_j` to
:math:`\bm{A}` keeps the photoionization it produces but places every
photoelectron at its point of emission, which removes the upstream seeding
that makes photoionization a feedback loop; in a strongly non-uniform gap
the ionization zone can be as thin as :math:`1/\kappa_j` even when
:math:`\kappa_j d \gg 1`, and the inception field then changes by tens of
per cent.  The stiff :math:`e^{\pm\kappa_j x}` modes this adds to
:math:`\bm{\mathcal{A}}` are handled by how the inception condition is
evaluated (:ref:`Chap:Numerics:Riccati`), not by the propagator.
