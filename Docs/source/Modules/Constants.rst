.. _Chap:Constants:

Constants.py
============

``Constants.py`` centralises the physical constants (Boltzmann constant,
elementary charge, speed of light) used throughout the project, sourced
from ``scipy.constants`` (CODATA 2018).  Every other module imports from
here instead of hard-coding numeric values.

.. literalinclude:: ../../../Constants.py
   :language: python

API reference
-------------

.. automodule:: Constants
   :members:
   :undoc-members:
