.. _Chap:Citing:

How to cite
===========

If you publish results obtained with ``Incept1D`` — inception voltages,
inception curves, growth rates, or transport tables exported from it — you
must cite the paper that presents the model:

   R. Marskar and C. Franck, *Role of negative-ion kinetics for electrical
   breakdown in air gaps*, J. Phys. D: Appl. Phys. (submitted, 2026).

.. code-block:: bibtex

   @article{Incept1D,
     author  = {Marskar, Robert and Franck, Christian},
     title   = {Role of negative-ion kinetics for electrical breakdown in air gaps},
     journal = {Journal of Physics D: Applied Physics},
     year    = {2026},
     note    = {Submitted}
   }

The same reference is in ``CITATION.cff`` at the root of the repository, so
GitHub's *Cite this repository* button gives it in APA and BibTeX form.

Cite the software as well, so that a reader can find the code.  Every
release is archived on Zenodo.
`doi:10.5281/zenodo.21916821 <https://doi.org/10.5281/zenodo.21916821>`_ is
the *concept* DOI: it stands for all versions and always resolves to the
latest.  Each release also has a DOI of its own, listed on that page, and
that is the one that identifies the version you used.  If you ran a commit
between releases, cite the concept DOI and give the git commit recorded in
the header of every output file (:ref:`Chap:Obtaining`), which identifies
the version precisely.

The swarm data shipped with the air mechanisms carry their own citation
requirements; cite the database named in the header of the file you used.
