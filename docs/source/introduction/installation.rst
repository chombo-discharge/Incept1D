.. _Chap:Installation:

How to install
==============

Incept1D is a pure-Python package with a ``pyproject.toml``; there is
nothing to compile.  An *editable* install (``-e``) is recommended: the
``incept1d`` command then always runs the sources in your checkout, so
pulling new commits or editing the code needs no reinstall.  The runtime
dependencies (:ref:`Chap:Prerequisites`) are installed automatically.

Using a virtual environment (recommended)
-----------------------------------------

.. code-block:: bash

   cd Incept1D
   python3 -m venv .venv
   source .venv/bin/activate          # on Windows: .venv\Scripts\activate
   pip install -e .                   # runtime dependencies + the incept1d command
   pip install -e .[docs]             # add Sphinx, only to build the docs
   pip install -e .[dev]              # add pytest, pre-commit, black, flake8

Into the user site
------------------

Without a virtual environment, ``pip install --user -e .`` installs the
command as ``~/.local/bin/incept1d`` (make sure that directory is on your
``PATH``).  On Debian/Ubuntu ≥ 23.04 pip refuses to install into a
system-managed interpreter (PEP 668); either use the virtual environment
above or add ``--break-system-packages`` — the editable install only adds a
path file to ``~/.local`` and can be undone with ``pip uninstall incept1d``.

Using system packages
---------------------

On Debian/Ubuntu the runtime dependencies are also available as system
packages, which pip will then leave alone:

.. code-block:: bash

   sudo apt install python3-numpy python3-scipy python3-matplotlib
   pip install --user --no-deps -e .

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

   incept1d --version
   incept1d pdiv mechanisms/air/air_pancheshnyi.py --pd-num 10 --no-plot

should print a small table of inception fields and voltages for a uniform
gap in dry air at 1 bar and exit without errors.  If the mechanism fails to
load with an ``AttributeError`` mentioning ``numpy``, check the version
mismatch note above.

Running the commands
--------------------

Every tool is a subcommand of ``incept1d`` (``incept1d --help`` lists them).
Mechanism and configuration files are given by path, so the commands can be
run from any directory; the examples in this documentation assume the
repository root:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py [CONFIG.json ...] [options]

Every subcommand accepts ``--help``.  See :ref:`Chap:QuickStart` for a
first walk-through and :ref:`Chap:ModulesOverview` for the full description
of each command.
