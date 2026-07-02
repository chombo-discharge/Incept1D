Introduction
============

Motivation
----------

Paschen's law states that the DC breakdown voltage of a uniform-field gap is
a universal function of the reduced gap distance :math:`pd`, where :math:`p`
is the gas pressure and :math:`d` is the electrode separation. It follows
from Townsend's ionization theory: avalanche growth is parameterized by the
first Townsend coefficient :math:`\alpha(E/p)` together with a cathode
secondary-electron-emission (SEE) efficiency :math:`\gamma`, giving the
classical breakdown criterion :math:`\alpha d = \ln(1 + \gamma^{-1})`.

Paschen's law is known to be asymptotically inaccurate beyond
:math:`pd \gtrsim 1-10~\mathrm{bar\,mm}`. In electronegative gases such as
air, part of that inaccuracy comes from a mechanism that the classical theory
omits entirely: **negative-ion transit and detachment**. Electron attachment
does not necessarily remove electrons permanently — if a negative ion drifts
across the gap and only *later* undergoes collisional or associative
detachment, it hands its electron back to the swarm mid-gap, where it can
keep contributing to avalanche growth. Whether this matters depends on how
the negative-ion transit time across the gap compares with the detachment
time, which is exactly the kind of question a purely local (`E/N`-only)
ionization coefficient cannot answer.

What this program computes
---------------------------

Incept1D implements a one-dimensional, multi-species drift-reaction model of
a discharge gap (electrons, positive ions, several negative-ion species, and
a two-stream photoionization field), and reduces the question "does a
self-sustained discharge start in this gap?" to a boundary-value / eigenvalue
problem that can be solved numerically for arbitrary pressure, gap length,
gas mixture, plasma chemistry, and electrode geometry (uniform field,
sphere-plane, or sphere-sphere).

Concretely, the code lets you:

- compute generalized Paschen curves (breakdown voltage / reduced field vs.
  :math:`pd`) that self-consistently include ion drift, attachment,
  detachment, ion-ion conversion, and photoionization/photoemission
  (:mod:`Inception`);
- inspect the *local* eigenvalues of the reaction-transport matrix
  :math:`R V^{-1}` as a function of :math:`E/N`, i.e. the purely local
  (non-transit) growth/decay modes of the plasma chemistry
  (:mod:`Eigenvalues`);
- compare the full boundary-value inception criterion against the classical
  ionization-integral criterion :math:`\int \max(\alpha - \eta, 0)\,dx`
  (:mod:`IonizationIntegral`);
- compute the temporal growth rate :math:`\lambda` of the discharge once the
  applied voltage exceeds the inception voltage (:mod:`Lambda`);
- export transport-coefficient and rate-coefficient tables from a mechanism
  file for use in external 3-D plasma solvers such as ``chombo-discharge``
  (:mod:`CreateChomboDischargeData`).

The full theoretical derivation is given in the accompanying manuscript,
``concepts.tex`` (section "Methods" → "Theoretical model"); :doc:`model`
summarizes the parts of that derivation needed to read the code, and
:doc:`modules/index` maps each equation onto the function that implements
it.

How the pieces fit together
----------------------------

The program has no package/build system: every script is run directly from
the repository root, e.g.

.. code-block:: console

   python Inception.py Air/Air_Pancheshnyi.py --p 1.0 --pd-min 1e-2 --pd-max 1e3

A *mechanism file* (e.g. ``Air/Air_Pancheshnyi.py``) declares the plasma
chemistry — species, reactions, transport coefficients, SEE efficiencies,
and photoionization data — behind a small fixed interface. ``Inception.py``
loads a mechanism file via :func:`Inception.load_mechanism` and exposes it as
a :class:`Inception.Mechanism` object; every other script
(:mod:`Eigenvalues`, :mod:`IonizationIntegral`, :mod:`Lambda`,
:mod:`CreateChomboDischargeData`) consumes that same object. This means a
new gas mixture, a new cross-section set, or a new set of reaction rates only
ever needs to be added in one place — a new mechanism file — and every tool
in the repository can immediately use it.
