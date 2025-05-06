.. include:: ../additional/links.rst
.. _isotropic-label:

Isotropic system
================

Here, we illustrate how the ``NMRDfromMD`` package can be applied to a
simple MD simulation. The NMR relaxation times :math:`T_1` and :math:`T_2`
are measured from a bulk polymer-water mixture using ``NMRDfromMD``.
To follow the tutorial, |MDAnalysis|, |NumPy|, and
|Matplotlib| must be installed.

MD system
---------

.. image:: ../figures/tutorials/isotropic-systems/snapshot-dark.png
    :class: only-dark
    :alt: PEG-water mixture simulated with LAMMPS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

.. image:: ../figures/tutorials/isotropic-systems/snapshot-light.png
    :class: only-light
    :alt: PEG-water mixture simulated with LAMMPS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

The system consists of a bulk mixture of 320 :math:`\text{H}_2\text{O}` water
molecules and 32 :math:`\text{PEG 300}` polymer molecules. The :math:`\text{TIP4P}-\epsilon`
is used for the water :cite:`fuentes-azcatlNonPolarizableForceField2014`.
:math:`\text{PEG 300}` refers to polyethylene glycol chains with a molar mass of
:math:`300~\text{g/mol}`. The trajectory was recorded during a
:math:`10~\text{ns}` production run performed using the open-source code
LAMMPS in the :math:`NpT` ensemble with a timestep of :math:`1~\text{fs}`.
The temperature was set to :math:`T = 300~\text{K}` and the pressure to
:math:`p = 1~\text{atm}`. Atomic positions were saved in the **prod.xtc** file
every :math:`1~\text{ps}`.

.. admonition:: Note
    :class: non-title-info

    If you'd like to learn LAMMPS and generate your own trajectories, beginner
    |lammps-tutorials| are available here :cite:`gravelleSetTutorialsLAMMPS2025`.

File preparation
----------------

To access the LAMMPS input files and pre-computed trajectory
data, download the |zip-peg-water-mixture| archive, or clone the
|dataset-peg-water-mixture| repository using:

.. code-block:: bash

    git clone https://github.com/NMRDfromMD/dataset-peg-water-mixture.git

The necessary trajectory files for this tutorial are located in the ``data/``
directory.

Create a MDAnalysis universe
----------------------------

Open a new Python script or Notebook, and define the path to the data files:

.. code-block:: python

	datapath = "mypath/polymer-in-water/data/"

Then, import ``NumPy``, ``MDAnalysis``, and the ``NMRD``
module of ``NMRDfromMD``:

.. code-block:: python

	import numpy as np
	import MDAnalysis as mda
	from nmrdfrommd import NMRD

From the trajectory files, create a ``universe`` by loading the
configuration file and trajectory:

.. code-block:: python

    u = mda.Universe(datapath+"production.data",
                     datapath+"production.xtc")

.. admonition:: Note
    :class: non-title-info
        
    The MDAnalysis ``universe``, ``u``, contains both the topology
    (atom types, masses, etc.) and the trajectory (atom positions
    at each frame). These informations are used by ``NMRDfromMD``
    for the calculation of NMR properties.

Let us print some basic information from the ``universe``, such as the number
of molecules (water and PEG):

.. code-block:: python

    n_TOT = u.atoms.n_residues
    n_H2O = u.select_atoms("type 6 7").n_residues
    n_PEG = u.select_atoms("type 1 2 3 4 5").n_residues
    print(f"The total number of molecules is {n_TOT} ({n_H2O} H2O, {n_PEG} PEG)")

Executing the script using Python will return:

.. code-block:: bw

    The total number of molecules is 450 (420 H2O, 30 PEG)

Let us also print information concerning the trajectory,
namely the timestep, ``timestep`` and
the total duration of the simulation, ``total_time``:

.. code-block:: python

    timestep = np.int32(u.trajectory.dt)
    total_time = np.int32(u.trajectory.totaltime)

    print(f"The timestep is {timestep} ps")
    print(f"The total simulation time is {total_time//1000} ns")

Executing the script using Python will return:

.. code-block:: bw

    The timestep is 2 ps
    The total simulation time is 10 ns

.. admonition:: Note
    :class: non-title-info

    In the context of ``MDAnalysis``, the ``timestep`` refers to the duration between
    two recorded frames, which is different from the actual timestep of
    :math:`1\,\text{fs}` used in the LAMMPS molecular dynamics simulation.

Launch the H-NMR analysis
-------------------------

Let us create three atom groups: the hydrogen atoms of the PEG, the hydrogen
atoms of the water, and all hydrogen atoms:

.. code-block:: python

    H_PEG = u.select_atoms("type 3 5")
    H_H2O = u.select_atoms("type 7")
    H_ALL = H_PEG + H_H2O

Then, let us first run NMRDfromMD for all hydrogen atoms:

.. code-block:: python

    nmr_all = NMRD(
        u=u,
        atom_group=H_ALL,
        neighbor_group = H_ALL,
        number_i=40)
    nmr_all.run_analysis()

With ``number_i = 40``, only 40 randomly selected atoms from ``H_ALL`` are
used in the calculation. Increase this number for better statistical resolution,
or set ``number_i = 0`` to include all atoms in the group. Here, ``H_ALL``
is specified as both the ``atom_group`` and ``neighbor_group``.

Run the script using Python.

Extract the NMR spectra
-----------------------

Let us access the calculated value of the NMR relaxation time :math:`T_1`
in :math:`f \to 0` by adding the following lines to the Python script:

.. code-block:: python

    T1 = np.round(nmr_all.T1, 2)

    print(f"The NMR relaxation time is T1 = {T1} s")

which should return:

.. code-block:: bw

    The NMR relaxation time is T1 = 1.59 s

The exact value you obtain will be different, as it depends on which hydrogen
atoms were randomly selected for the calculation. With the small value
``number_i = 40``, the noise is important. You can increase that number
for more precise result, but this will obviously increase the computation time.

The full :math:`T_1` and :math:`T_2` spectra can be extracted as
``1/nmr_ALL.R1`` and ``1/nmr_ALL.R2``, respectively. The corresponding
frequencies are stored in ``nmr_ALL.f``.

.. code-block:: python

    R1_spectrum = nmr_ALL.R1
    R2_spectrum = nmr_ALL.R2
    T1_spectrum = 1 / R1_spectrum
    T2_spectrum = 1 / R2_spectrum
    f = nmr_ALL.f

The spectra :math:`T_1` and :math:`T_2` can then be plotted as a function of
:math:`f` using ``pyplot``:

.. code-block:: python

    from matplotlib import pyplot as plt
    plt.loglog(f, T1_spectrum, 'o', label='T1')
    plt.loglog(f, T2_spectrum, 's', label='T2')
    plt.xlabel("f (MHz)")
    plt.ylabel("T1, T2 (s)")
    plt.legend()
    plt.show()

.. image:: ../figures/tutorials/isotropic-systems/T1-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/tutorials/isotropic-systems/T1-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: NMR relaxation times :math:`T_1` (disks) and :math:`T_2` (squares)
    as a function of the frequency :math:`f` for the
    :math:`\text{PEG-H}_2\text{O}` bulk mixture.

Calculate the intra-molecular NMR relaxation
--------------------------------------------

Let us measure the intra-molecular contribution to the NMR relaxation time.
To simplify the analysis, we will differentiate between PEG and
:math:`\text{H}_2\text{O}` molecules and perform two separate analyses.

.. code-block:: python

    nmr_H2O_intra = nmrmd.NMR(u, atom_group = H_H2O, type_analysis="intra_molecular", number_i=40)
    nmr_PEG_intra = nmrmd.NMR(u, atom_group = H_PEG, type_analysis="intra_molecular", number_i=40)

The correlation function :math:`G_{ij}` can be accessed from
``nmr_H2O_intra.gij[0]``, and the corresponding time values from
``nmr_H2O_intra.t``:

.. code-block:: python

    t = nmr_PEG_intra.t
    G_intra_H2O = nmr_H2O_intra.gij[0]
    G_intra_PEG = nmr_PEG_intra.gij[0]

.. image:: ../figures/tutorials/isotropic-systems/Gintra-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water-PEG

.. image:: ../figures/tutorials/isotropic-systems/Gintra-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water-PEG

.. container:: figurelegend

    Figure: Intra-molecular correlation function :math:`G_\text{R}` for both
    PEG (squares) and :math:`\text{H}_2\text{O}` (disks).

From the correlation functions, one can obtain the typical rotational
time of the molecules:

.. code-block:: python

    tau_rot_H2O = np.round(np.trapz(G_intra_H2O, t)/G_intra_H2O[0],2)
    tau_rot_PEG = np.round(np.trapz(G_intra_PEG, t)/G_intra_PEG[0],2)
    print(f"The rotational time of H2O is = {tau_rot_H2O} ps")
    print(f"The rotational time of PEG is = {tau_rot_PEG} ps")

This should return (again, the exact values will likely differ in your
case):

.. code-block:: bw

    >> The rotational time of H2O is = 6.35 ps
    >> The rotational time of PEG is = 8.34 ps