Eigenvalues
===========

``Eigenvalues.py`` computes the eigenvalues of the *local* transport matrix
:math:`\bm{A} = \bm{R}\bm{V}^{-1}` as a function of :math:`E/N`, at a single
spatial point — i.e. without any gap integration or boundary conditions.
This is the purely local counterpart of the full inception criterion in
:mod:`Inception`: a positive real eigenvalue means that, *locally*, net
ionization exceeds attachment/loss and the charge-carrier flux would grow
exponentially in an infinite, spatially uniform medium. As shown in
:doc:`../model`, the largest such eigenvalue reduces to :math:`\lambda_+`
of the closed-form three-species model when the mechanism only has
ionization/attachment/detachment, and it is the "apparent effective
ionization coefficient" referenced throughout the manuscript.

Because eigenvalues have no inherent ordering, :func:`_track_step` uses the
Hungarian algorithm to match each eigenvalue at one :math:`E/N` sample to
its continuation at the next sample (minimizing total squared distance in
the complex plane), so that plotted eigenvalue "tracks" don't jump between
unrelated modes as :math:`E/N` varies.

.. literalinclude:: ../../Eigenvalues.py
   :language: python
   :pyobject: compute_eigenvalues

Two CLI modes are available: the default mode plots every eigenvalue track
vs. :math:`E/N` at a fixed pressure; ``--pressure-scan`` instead plots one
selected eigenvalue track (default: the leading/most-growing mode) vs.
:math:`E/N` for several pressures, which is how the manuscript's discussion
of pressure scaling of the detachment process was produced.

Command-line usage
--------------------

.. code-block:: console

   python Eigenvalues.py <mechanism.py> [CONFIG.json ...] \
       [--p P] [--T T] [--EN-lo LO] [--EN-hi HI] [--EN-num N] \
       [--pressure-scan] [--p-min MIN] [--p-max MAX] [--p-num N] \
       [--eig-index IDX] [--write-to-file FILE]

API reference
--------------

.. automodule:: Eigenvalues
   :members:
   :undoc-members:
