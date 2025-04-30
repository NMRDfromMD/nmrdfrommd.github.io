.. include:: ../additional/links.rst
.. _anisotropic-label:

Anisotropic systems
===================

In this tutorial, the NMR relaxation rate :math:`R_1` is measured from water
confined in a nanoslit of silica.

.. admonition:: Note
    :class: non-title-info
    
    If you are new to NMRDfromMD, it is recommended that you follow this simpler
    tutorial on :ref:`isotropic-label` first.

MD system
---------

.. container:: hatnote

   Measuring the NMR relaxation time of nanoconfined water

.. image:: ../figures/tutorials/anisotropic-systems/snapshot-dark.png
    :class: only-dark
    :alt: Water confined in silica slit with GROMACS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

.. image:: ../figures/tutorials/anisotropic-systems/snapshot-light.png
    :class: only-light
    :alt: Water confined in silica slit with GROMACS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

The system is composed of 602 :math:`\text{TIP4P}-\epsilon` water molecules
confined in a slit silica nanopore. The trajectory was recorded during a
:math:`10\,\text{ns}` production run performed using the open-source code GROMACS,
in the anisotropic :math:`NP_zT` ensemble, with a timestep of :math:`1\,\text{fs}`.
To balance the surface charge, 20 sodium ions are present in the slit.
The imposed temperature was :math:`T = 300\,\text{K}`, and the pressure
was :math:`p = 1\,\text{bar}`. Atomic positions were saved in the **prod.xtc**
file every :math:`2\,\text{ps}`.

.. admonition:: Note
    :class: non-title-info

    If you'd like to learn GROMACS and generate your own trajectories, beginner
    |gromacs-tutorials| are available here.

File preparation
----------------

To access all trajectory and input files, download the
``water-in-silica`` repository from GitHub, or simply clone it
to your computer using:

.. code-block:: bash

    git clone https://github.com/simongravelle/water-in-silica.git

The dataset required to follow this tutorial is located in
``raw-data/N50/``.

Create a MDAnalysis universe
----------------------------

Open a new Python script or a new notebook, and define
the path to the data files:

.. code-block:: python

    datapath = "mypath/water-in-silica/raw-data/N50/"

Then, import NumPy, MDAnalysis, and NMRforMD:

.. code-block:: python

    import numpy as np
    import MDAnalysis as mda
    import nmrformd as nmrmd

From the trajectory files, create a MDAnalysis universe by
loading the configuration file and the trajectory:

.. code-block:: python

    u = mda.Universe(datapath + "prod.tpr", datapath + "prod.xtc")

Next, extract some basic information from the universe, such as
the number of molecules, the timestep, and the total duration:

.. code-block:: python

    n_molecules = u.atoms.n_residues
    print(f"The number of molecules is {n_molecules}")
    timestep = np.int32(u.trajectory.dt)
    print(f"The timestep is {timestep} ps")
    total_time = np.int32(u.trajectory.totaltime)
    print(f"The total simulation time is {total_time} ps")

Running this will return:

.. code-block:: bw

    >> The number of molecules is 623
    >> The timestep is 2 ps
    >> The total simulation time is 10000 ps

Launch the NMR analysis
-----------------------

Let us create three atom groups corresponding to: the hydrogen atoms of the silica,
the hydrogen atoms of the water, and all the hydrogen atoms combined:

.. code-block:: python

    H_H2O = u.select_atoms("name HW1 HW2")
    H_SIL = u.select_atoms("name H")
    H_ALL = H_H2O + H_SIL

Then, let us run three separate NMR analyses: one for the water-silica interaction only,
one for the intra-molecular interaction of water, and one for the inter-molecular
interaction of water:

.. code-block:: python

    nmr_H2O_SIL = nmrmd.NMR(u, atom_group = H_H2O,
                        neighbor_group = H_SIL, number_i=40, isotropic=False)
    nmr_H2O_INTRA = nmrmd.NMR(u, atom_group = H_H2O, neighbor_group = H_H2O, number_i=40,
                            type_analysis = 'intra_molecular', isotropic=False)
    nmr_H2O_INTER = nmrmd.NMR(u, atom_group = H_H2O, neighbor_group = H_H2O, number_i=40,
                            type_analysis = 'inter_molecular', isotropic=False)

Note the use of ``isotropic=False``, which is required here because the system is anisotropic.

Extract the NMR spectra
-----------------------

Let us access the NMR relaxation rate :math:`R_1`:

.. code-block:: python

    R1_spectrum_H2O_SIL = nmr_H2O_SIL.R1
    R1_spectrum_H2O_INTRA = nmr_H2O_INTRA.R1
    R1_spectrum_H2O_INTER = nmr_H2O_INTER.R1
    f = nmr_H2O_SIL.f

.. image:: ../figures/tutorials/anisotropic-systems/spectra-dark.png
    :class: only-dark
    :alt: NMR results obtained from the GROMACS simulation of water in silica

.. image:: ../figures/tutorials/anisotropic-systems/spectra-light.png
    :class: only-light
    :alt: NMR results obtained from the GROMACS simulation of water in silica

.. container:: figurelegend

    Figure: NMR relaxation rates :math:`R_1` for the water confined in
    a silica slit.

Note that the :math:`\text{H}_2\text{O}-\text{silica}` contribution is much
smaller than the intra- and intermolecular contributions from the water. This
is due to the relatively small number of hydrogen atoms from the silica (92),
compared to the 1204 hydrogen atoms from the water.
