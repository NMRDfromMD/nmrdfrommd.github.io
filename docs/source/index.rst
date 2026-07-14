.. include:: additional/links.rst

NMRDfromMD
==========

NMRDfromMD is a Python toolkit for computing nuclear magnetic resonance (NMR)
relaxation properties directly from molecular dynamics (MD) trajectories.
From atomistic simulations, it evaluates the dipolar interactions between
nuclear spins to calculate the longitudinal (:math:`T_1`) and transverse
(:math:`T_2`) relaxation times. NMR relaxation is sensitive to molecular
motion over a broad range of timescales, making it a powerful probe of
translational and rotational dynamics in liquids, confined fluids, polymers,
and biological systems. By extracting relaxation times from MD simulations,
NMRDfromMD enables direct comparison with experimental NMR measurements,
allowing users to validate simulation models and identify the molecular
mechanisms underlying relaxation. It can also be used to predict
and interpret NMR relaxation behavior when experimental data are unavailable.

Compatible Simulation Packages
-------------------------------

NMRDfromMD accepts any trajectory format supported by MDAnalysis,
covering virtually all major MD simulation packages including
|LAMMPS|, |GROMACS|, NAMD, AMBER, CHARMM, and many others.
For a full list of supported formats, see the |MDAnalysis| documentation.

This package builds on the now discontinued |NMRforMD|.

.. image:: ../../avatars/avatars.png
    :class: only-dark
    :alt: molecular dynamics systems used in these examples 

.. image:: ../../avatars/avatars.png
    :class: only-light
    :alt: molecular dynamics systems used in these examples

.. container:: figurelegend

    Figure: Examples of systems that can be analyzed with NMRDfromMD, spanning
    simple bulk liquids (water, left), idealized Lennard-Jones fluids (center),
    and biologically relevant systems such as a lysozyme protein and its
    hydration shell (right).

Installation
------------

To install the latest development version of NMRDfromMD, clone the repository,
|NMRDfromMD-code|, from GitHub, and use ``pip`` from the main directory:

.. code-block:: bash

    git clone https://github.com/NMRDfromMD/nmrdfrommd.git

    cd nmrdfrommd/

    pip install .


To run the test suite, install the testing dependencies:

.. code-block:: bash

    pip install pytest coverage

Then run from the tests folder:

.. code-block:: bash

    pytest

Datasets
--------

Molecular dynamics datasets are available on GitHub: a |polymer in water|
system generated using LAMMPS, and a |water confined in silica| system
generated using GROMACS. These datasets can be downloaded to follow the
tutorials or simply to test NMRDfromMD.

.. toctree::
   :maxdepth: 2
   :caption: Theory
   :hidden:

   theory/context
   theory/mathematical-framework
   applications/lennard-jones-fluids

.. toctree::
   :maxdepth: 2
   :caption: Getting started
   :hidden:

   applications/tutorial
   theory/best-practice

.. toctree::
   :maxdepth: 2
   :caption: Applications
   :hidden:

   applications/lysozyme-in-water
   applications/anisotropic-system

.. toctree::
   :maxdepth: 2
   :caption: Additional
   :hidden:

   theory/microscopic-origin
   additional/scope
   additional/simulation-methods
   additional/bibliography
   additional/acknowledgments
