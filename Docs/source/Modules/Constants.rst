.. _Chap:Constants:

Physical constants — ``incept1d.constants``
===========================================

:mod:`incept1d.constants` centralises the physical constants (Boltzmann constant,
elementary charge, speed of light) used throughout the project, sourced
from ``scipy.constants`` (CODATA 2018).  Every other module imports from
here instead of hard-coding numeric values.

.. literalinclude:: ../../../src/incept1d/constants.py
   :language: python

API reference
-------------

.. automodule:: incept1d.constants
   :members:
   :undoc-members:
