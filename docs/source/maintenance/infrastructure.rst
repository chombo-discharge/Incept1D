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

The figures in :ref:`Chap:Examples:IEC60052`, :ref:`Chap:Examples:Electra`
and :ref:`Chap:Photoionization` are **built, not shipped**.  Each one is
the product of three committed ingredients:

1. **Reference data** — the IEC 60052 tables and the ELECTRA compilation in
   ``examples/``, with a header giving their provenance.
2. **The calculation** — ``examples/*/run.sh`` (which call
   ``incept1d pdiv``) and ``mechanisms/air/zheleznyak.py``, producing ``.dat`` files
   with a metadata header that records the git commit and the full
   command line.
3. **The plot** — a pgfplots source in ``docs/figures/`` that reads the
   ``.dat`` files by column and is compiled with ``pdflatex``.

``docs/figures/Makefile`` chains them: it copies any reference data and
runs the calculations into ``docs/build/figures/<Figure>/``, compiles the
``.tex`` there, and places a PDF (for the LaTeX builder) and a 150 dpi PNG
(for the HTML builder) in ``docs/source/figures/``, which is ignored by
git.  The Sphinx sources reference them as ``figures/<Figure>.*`` so each
builder picks its format.

.. code-block:: bash

   cd Docs
   make figures                 # everything
   make -C figures Electra      # one figure
   PD_NUM=30 make figures       # quicker, coarser curves for a local check

The computed curves are only re-run when something that can change them
is newer: the solver modules, anything under ``mechanisms/air/``, or the run scripts.
The full IEC set takes tens of minutes; in CI the computed ``.dat`` files
are cached under a key derived from those same inputs, so a
documentation-only change rebuilds in minutes while a physics change
recomputes the curves.  A figure whose plot overflows its page (pgfplots
silently spills onto a second page) fails the build.

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
