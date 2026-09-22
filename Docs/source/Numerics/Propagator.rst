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
The full augmented matrix is small (:math:`N_s + 2N_\gamma`, i.e. 6–18 for
the air schemes), so a single ``expm`` costs microseconds.

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

The initial :math:`N_{\min}` is additionally raised by
:func:`FieldDistributions._compute_n_steps` for sphere geometries so that
the field does not change by more than a set fraction over the first step
near the sphere, where :math:`f(\xi)` is steepest.

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

Local treatment of strongly absorbed photons
--------------------------------------------

Photon groups whose optical depth over the gap is large,
:math:`\kappa_j d > 12`, are re-absorbed within a small fraction of any
integration step.  Propagating them explicitly adds stiff decay modes to
:math:`\bm{\mathcal{A}}` without changing the physics, so
:func:`incept1d.solver._build_A_aug` removes such groups from the augmented
system and adds their steady-state contribution
:math:`2\vec{b}_j\vec{c}_j^\intercal/\kappa_j` directly to :math:`\bm{A}`
(the factor 2 accounts for both streams).  The number of explicitly
propagated groups :math:`N_\gamma^\mathrm{eff}` is decided once per
determinant evaluation and used consistently in :math:`\bm{Q}`.
