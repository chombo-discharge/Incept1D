# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Incept1D — discharge inception in a 1-D drift-reaction model.

Library layout
--------------
``constants``    physical constants
``reactions``    reaction-string parser → reaction-rate matrix R
``fields``       gap geometry / normalised field profiles f(ξ)
``mechanism``    loading mechanism files and JSON configurations
``solver``       augmented ODE, propagators, boundary determinant det Q(λ)
``inception``    roots of det Q in E/N and inception-curve branch tracking
``eigenvalues``  eigenvalues of the local transport matrix R V⁻¹
``ionization``   ionization integrals across the gap
``growth``       temporal growth rate λ above inception
``chombo``       transport tables for the chombo-discharge 3-D solver
``output``       shared result-file metadata header
``cli``          the ``incept1d`` command-line interface
"""

__version__ = "0.1.0"
