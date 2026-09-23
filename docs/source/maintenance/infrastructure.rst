.. _Chap:Infrastructure:

Hosting, CI and documentation builds
====================================

.. contents:: On this page
   :local:
   :depth: 1

Hosting
-------

The code is hosted on GitHub at |repo| under the ``chombo-discharge``
organisation; issues and pull requests are handled there, and the
documentation is deployed to
`chombo-discharge.github.io/Incept1D <https://chombo-discharge.github.io/Incept1D/>`_.
Contributors work on forks and open pull requests against ``main``.
There is no package on PyPI; the repository is the distribution.

Continuous integration
----------------------

GitHub Actions workflows under ``.github/workflows/`` run on every push
and pull request:

* **docs** (``docs.yml``) — builds the figures and the HTML documentation
  with warnings treated as errors, uploads the result as a build artifact,
  and — on ``main`` — deploys that artifact to GitHub Pages.  Nothing
  pre-rendered is committed: the deployed site is entirely a product of
  the repository at that commit.
* **Lint** and **tests** (planned, see :ref:`Chap:TestSuite`) —
  ``pre-commit run --all-files`` and ``python3 -m pytest tests/``.

The Sphinx build is also part of the pre-commit configuration (as a
``dummy`` build that runs the whole pipeline without writing output), so
broken cross-references or autodoc import failures are caught before a
commit is made.

How the figures are built
-------------------------

The figures in :ref:`Chap:Examples:Paschen` and :ref:`Chap:Photoionization`
are **built, not shipped**.  Each is the product of two committed
ingredients:

1. **The calculation** — ``examples/*/run.sh`` (which call ``incept1d pdiv``)
   and ``mechanisms/air/zheleznyak.py``, producing ``.dat`` files with a
   metadata header that records the git commit and the full command line.
2. **The plot** — a pgfplots source in ``docs/figures/`` that reads the
   ``.dat`` files by column and is compiled with ``pdflatex``.

``docs/figures/Makefile`` chains them: it runs the calculations into
``docs/build/figures/<figure>/``, compiles the ``.tex`` there, and places a
PDF (for the LaTeX builder) and a 150 dpi PNG (for the HTML builder) in
``docs/source/figures/``, which is ignored by git.  The Sphinx sources
reference them as ``figures/<figure>.*`` so each builder picks its format.

.. code-block:: bash

   cd docs
   make figures                   # the default figures (seconds)
   make -C figures paschen        # one figure
   PD_NUM=30 make figures         # quicker, coarser curves for a local check

Default and opt-in figures
..........................

``make figures`` builds only ``paschen`` and ``zheleznyakfit``.  Both are
self-contained and take seconds.

Figures that compare against **published reference data** are **opt-in**
(``OPTIONAL_FIGURES`` in the Makefile), for two reasons.  Such data is
generally copyrighted and therefore not part of this repository, so without it
the figure would show the computed curves with nothing to compare against; and
these are the expensive figures, together the large majority of the total
figure cost.  A contributor who holds the data supplies it as the relevant
example page describes, and builds the figure by name:

.. code-block:: bash

   make -C docs/figures <figure>

The computed curves are only re-run when something that can change them is
newer: the solver modules, anything under ``mechanisms/``, or the run
scripts.  In CI the ``.dat`` files are cached under a key derived from those
same inputs.  A figure whose plot overflows its page (pgfplots silently spills
onto a second page) fails the build.

Requirements: ``pdflatex`` with ``pgfplots`` and ``siunitx`` (TeX Live
packages ``texlive-latex-extra``, ``texlive-pictures``,
``texlive-science``) and ``pdftoppm`` (``poppler-utils``).

Building the documentation locally
----------------------------------

The sources are in ``docs/source``; build products go to
``docs/build/<builder>``.

.. code-block:: bash

   cd Docs
   make html          # figures + HTML -> build/html/index.html
   make latexpdf      # figures + PDF  -> build/latex/Incept1D.pdf (needs xelatex)
   make dummy         # full Sphinx check, no output (what pre-commit runs)
   make clean         # removes build/ and the generated figures

``make`` invokes ``python3 -m sphinx`` rather than the bare
``sphinx-build`` executable, so that the interpreter that has NumPy, SciPy
and matplotlib installed is the one used for autodoc.  The equivalent
direct command from the repository root is

.. code-block:: bash

   python3 -m sphinx -W --keep-going -b html docs/source docs/build/html

Autodoc imports the project modules, so the runtime dependencies must be
installed in the same environment as Sphinx.  ``conf.py`` forces the
``Agg`` matplotlib backend and adds the repository root to ``sys.path``
for this purpose.

Conventions used in the sources
...............................

* One directory per section (``Introduction/``, ``Theory/``, ``Numerics/``,
  ``Modules/``, ``examples/``, ``Maintenance/``) plus ``figures/`` for the
  generated figures; the master ``index.rst``
  holds one hidden ``toctree`` per section so that the sidebar shows the
  sections as expandable headings.
* Every page starts with a label ``.. _Chap:Name:`` and is referenced with
  ``:ref:`Chap:Name```.
* Equations that are referenced get a ``:label:`` and are cited with
  ``:eq:``; figures and tables use ``:numref:``.
* Functions that directly implement an equation are shown with
  ``literalinclude`` (``:pyobject:``) rather than copied, so the docs
  cannot drift from the code.
* Literature is cited with ``[Key]_`` and listed in ``zzreferences.rst``.
  Sphinx fails the build on an unreferenced citation.
