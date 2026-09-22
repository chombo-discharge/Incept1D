.. _Chap:Installation:

How to install
==============

There is nothing to compile or install: the scripts are run directly from
the repository root with the system Python.  Only the dependencies listed in
:ref:`Chap:Prerequisites` must be available.

Using a virtual environment (recommended)
-----------------------------------------

.. code-block:: bash

   cd Incept1D
   python3 -m venv .venv
   source .venv/bin/activate          # on Windows: .venv\Scripts\activate
   pip install numpy scipy matplotlib
   pip install -r Docs/requirements.txt   # only if you want to build the docs
   pip install pre-commit pytest          # only for development

Using system packages
---------------------

On Debian/Ubuntu the runtime dependencies are available as

.. code-block:: bash

   sudo apt install python3-numpy python3-scipy python3-matplotlib

.. warning::

   Mixing ``apt``-installed and ``pip``-installed NumPy/matplotlib on the
   same interpreter is a common source of ``numpy.core.multiarray failed to
   import`` errors: a ``pip``-installed NumPy 2.x in ``/usr/local`` shadows
   the ``apt`` NumPy 1.x that the ``apt`` matplotlib was compiled against.
   If you see that error, either use a virtual environment (above) or
   upgrade matplotlib in the same layer as NumPy
   (``pip install --upgrade matplotlib contourpy kiwisolver``).

Verifying the installation
--------------------------

From the repository root:

.. code-block:: bash

   python3 Inception.py Air/Air_Pancheshnyi.py --pd-num 10 --no-plot

should print a small table of breakdown fields and voltages for a uniform
gap in dry air at 1 bar and exit without errors.  If the mechanism fails to
load with an ``AttributeError`` mentioning ``numpy``, check the version
mismatch note above.

Running the scripts
-------------------

All scripts are run from the repository root and take a mechanism file as
their first argument:

.. code-block:: bash

   python3 Inception.py Air/Air_Pancheshnyi.py [CONFIG.json ...] [options]

Every script accepts ``--help``.  See :ref:`Chap:QuickStart` for a first
walk-through and :ref:`Chap:ModulesOverview` for the full description of
each script.
