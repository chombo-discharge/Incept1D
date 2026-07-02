Constants
=========

``Constants.py`` centralizes the physical constants (Boltzmann constant,
elementary charge, speed of light) used throughout the project, sourced from
``scipy.constants`` (CODATA 2018). Every other module imports from here
instead of hard-coding numeric values, so a precision update only needs to
happen in one place.

.. literalinclude:: ../../Constants.py
   :language: python
   :lines: 7-9

API reference
--------------

.. automodule:: Constants
   :members:
   :undoc-members:
