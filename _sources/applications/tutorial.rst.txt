.. include:: ../additional/links.rst
.. _isotropic-label:

Tutorial
========

In this tutorial, the proton (:math:`^1\mathrm{H}`) NMR
relaxation properties of a polymer-water mixture is computed from a
molecular dynamics trajectory. More specifically,
NMRDfromMD is used to extract the frequency-dependent
relaxation rates [:math:`R_1(f)` and :math:`R_2(f)`], separate the
intra- and intermolecular contributions to the relaxation, and diferentiate 
the relative contribution from water and from the polymer to the relaxation.

To follow the tutorial, |MDAnalysis|, |NumPy|, and
|Matplotlib| must be installed alongside NMRDfromMD.

File preparation
----------------

The system consists of a bulk mixture containing 420 water
molecules and 30 polyethylene glycol (PEG 300) chains, simulated at :math:`T=300~\mathrm{K}`
and :math:`p=1~\mathrm{atm}` using LAMMPS. Details are given in :ref:`simulation-methods`.
To access the LAMMPS input files and pre-computed trajectory data, clone the
|dataset-peg-water-mixture| repository. The required trajectory files are
located in the ``data/`` directory.

.. admonition:: Important
    :class: warning

    The trajectory files are stored using Git Large File Storage (Git LFS).
    This means that after cloning the repository, you must download the
    actual trajectory data before using it.

    If Git LFS is not installed, install it first:

    .. code-block:: bash

        apt install git-lfs
        git lfs install

    Then retrieve the trajectory files:

    .. code-block:: bash

        git lfs pull

Alternatively, you can regenerate the trajectory by rerunning the LAMMPS
simulation scripts provided in the repository.

Import the simulation data into Python
--------------------------------------

Open a new Python script or Notebook, and define the path to the data files:

.. code-block:: python

	datapath = "mypath/dataset-peg-water-mixture/data/"

Then, import NumPy, MDAnalysis, and the NMRD module of NMRDfromMD:

.. code-block:: python

	import numpy as np
	import MDAnalysis as mda
	from nmrdfrommd import NMRD

From the trajectory files, create a MDAnalysis Universe by loading the
configuration file and trajectory:

.. code-block:: python

    u = mda.Universe(
        datapath+"production.data",
        datapath+"production.xtc"
    )

The Universe is a central object in MDAnalysis. It combines the system topology
(atom identities, masses, molecule definitions, etc.) with the time-dependent
atomic coordinates. NMRDfromMD uses this Universe to access both
the molecular structure and the atomic trajectories required to compute dipolar
correlation functions and NMR relaxation properties.

Before starting the analysis, it is useful to verify that the system has been
correctly loaded and to inspect its basic composition:

.. code-block:: python

    H2O = u.select_atoms("type 6 7")
    PEG = u.select_atoms("type 1 2 3 4 5")

    n_TOT = len(u.residues)
    n_H2O = len(H2O.residues)
    n_PEG = len(PEG.residues)

    print(
        f"System: {n_TOT} molecules"
        f"({n_H2O} H2O, {n_PEG} PEG)"
    )

Executing the script using Python will return:

.. code-block:: bw

    System: 450 molecules (420 H2O, 30 PEG)

This output confirms that the simulation contains the expected 450 molecules,
correctly partitioned into 420 water molecules and 30 PEG chains. This corresponds to
an ethylene oxide (EO) over water ratio of 0.46.

Let us also print information about the trajectory: ``frame_dt``
and the total duration of the simulation, ``total_time``:

.. code-block:: python

    frame_dt = np.int32(u.trajectory.dt)
    total_time = np.int32(u.trajectory.totaltime)

    print(f"Trajectory: {frame_dt} ps/frame, "
        f"{total_time//1000} ns total")

Executing the script using Python will return:

.. code-block:: bw

    Trajectory: 2 ps/frame, 10 ns total

The ``frame_dt`` refers to the time interval between two
stored frames in the trajectory file, not the integration timestep used
in the molecular dynamics simulation. In this tutorial, the trajectory was generated with LAMMPS using a
1 fs integration timestep, but atomic configurations were written to the
``production.xtc`` file every 2 ps. As a result, ``frame_dt`` corresponds to
the 2 ps sampling interval, which determines the
temporal resolution of all subsequent correlation functions and relaxation calculations.

Run the :math:`^1\text{H}` NMR relaxation analysis
--------------------------------------------------

First, three atom groups are created: one group containing the hydrogen atoms
belonging to PEG, another group containing hydrogen atoms belonging to water,
and the combined set containing all hydrogen atoms:

.. code-block:: python

    H_H2O = H2O.select_atoms("type 7")
    H_PEG = PEG.select_atoms("type 3 5")
    H_ALL = H_PEG + H_H2O

To extract the NMR relaxation properties from the entire system, NMRDfromMD is
first executed on all hydrogen atoms:

.. code-block:: python

    nmr_all = NMRD(
        u=u,
        atom_group=H_ALL,
        number_i=20
    )
    res_nmr = nmr_all.run_analysis()

On a standard laptop (Intel Core i9-12900H processor), this step
typically takes 1-2 minutes. The runtime depends mainly on three factors:
(i) the number of selected reference atoms (``number_i``),
(ii) the number of atoms in the system,
and (iii) the number of saved trajectory frames.

The parameter ``number_i`` controls how many reference hydrogen atoms are
randomly selected for the calculation. Computing the dipolar interaction
for every hydrogen atom can become computationally expensive in large
systems. Instead, ``NMRDfromMD`` samples a subset of reference atoms
while retaining all neighbouring atoms. This sampling reduces the computational
cost at the expense of increased statistical uncertainty.

All calculated values are stored within the ``res_nmr`` dictionary.
Let us first extract the NMR relaxation time :math:`T_1` in
the limit :math:`f \to 0`, add the following lines to the Python script:

.. code-block:: python

    T10 = res_nmr["T1"][0]
    T10_err = res_nmr["T1_err"][0]

    print(f"The NMR relaxation time is T1 = {T10:.2f} pm {T10_err:.2f} s")

The expected output is

.. code-block:: bw

    The NMR relaxation time is T1 = 0.64 pm 0.47 s

The error on :math:`T_1` is relatively large due to the small 
value of ``number_i`` used. Increasing
``number_i`` reduces this statistical uncertainty. Avalue of 200
for ``number_i`` returns:

.. code-block:: bw

    The NMR relaxation time is T1 = 0.64 pm 0.47 s

For reference, experimental measurements by Jora report the value for :math:`T_1`
for a water PEG 300 mixute with a similar ethylene oxide (EO) over water
ratio equal to :math:`T_1 = 0.7` s :cite:`joraDynamicalAspectsWaterpolyethylene2016`.

Alongside the relaxation time, :math:`T_1`, the longitudinal and transverse
relaxation rates,

.. math::

   R_1(f)=\frac{1}{T_1(f)}, \qquad
   R_2(f)=\frac{1}{T_2(f)},

are available for every frequency :math:`f` (in MHz) as ``res_nmr["R1"]``
and ``res_nmr["R2"]``. The corresponding frequencies are stored in
``res_nmr["f"]``.

.. code-block:: python

    R1_spectrum = res_nmr["R1"]
    R2_spectrum = res_nmr["R2"]
    f = res_nmr["f"]

The spectra :math:`R_1 (f)` and :math:`R_2 (f)` can then be plotted as a
function of :math:`f` using pyplot:

.. code-block:: python

    from matplotlib import pyplot as plt

    # Plot settings
    plt.figure(figsize=(8, 5))
    plt.loglog(
        f,
        R1_spectrum,
        "o",
        label="R1",
        markersize=5,
    )
    plt.loglog(
        f,
        R2_spectrum,
        's',
        label='R2',
        markersize=5
    )
    # Labels and Title
    plt.xlabel(
        "Frequency (MHz)",
        fontsize=12
    )
    plt.ylabel(
        "Relaxation Rates (s-1)",
        fontsize=12
    )
    # Grid and boundaries
    plt.grid(
        True, 
        which="both", 
        linestyle='--', 
        linewidth=0.7
    )
    plt.xlim(80, 1e5)
    plt.ylim(0.05, 2)
    # Legend
    plt.legend()
    plt.tight_layout()
    plt.show()

The resulting spectra should resemble :ref:`Fig. 1 <fig:nmr-tutorial1>`, panel A. For an
isotropic liquid, :math:`R_1(f)` and :math:`R_2(f)` are expected to
approach similar values in the low-frequency limit. In this regime, both relaxation
rates probe the low-frequency limit of the spectral density, which is
dominated by long-time molecular reorientations and translational diffusion.

Because only ``number_i = 20`` reference atoms are sampled here, the
spectra exhibit noticeable statistical noise. Repeating the calculation
with a larger value for ``number_i`` produces much smoother curves, as
shown :ref:`Fig. 1 <fig:nmr-tutorial1>`, panel B.

.. _fig:nmr-tutorial1:

.. container:: figure

    .. image:: tutorial/nmr-total-dm.png
        :class: only-dark
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. image:: tutorial/nmr-total.png
        :class: only-light
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. container:: figurelegend

        Figure 1: (A) :math:`^1\text{H}`-NMR relaxation
        rates :math:`R_1` and :math:`R_2` as a
        function of the frequency :math:`f` for the 
        :math:`\text{PEG-H}_2\text{O}` bulk mixture. Results are given for
        a small value of ``number_i``, :math:`n_i = 20`.
        (B) Same quantity as in panel A, but performed using all hydrogen atoms
        using :math:`n_i = 0`.

Separating intra and intermolecular
-----------------------------------

So far, the relaxation rates were calculated without distinguishing
between intra- and intermolecular interactions. One of the major
advantages of molecular dynamics simulations is that every dipolar
interaction can be classified according to whether it originates from two
nuclei within the same molecule or from two different molecules. This
separation is generally not accessible from experimental measurements
alone.

Let us extract the intramolecular contributions to the relaxation for 
both water and PEG, separately:

.. code-block:: python

    nmr_h2o_intra = NMRD(
        u=u,
        atom_group=H_H2O,
        type_analysis="intra_molecular",
        number_i=200)
    res_h2o_intra = nmr_h2o_intra.run_analysis()

    nmr_peg_intra = NMRD(
        u=u,
        atom_group=H_PEG,
        type_analysis="intra_molecular",
        number_i=200)
    res_peg_intra = nmr_peg_intra.run_analysis()

Similarly, we can also measure the intermolecular contributions:

.. code-block:: python

    nmr_h2o_inter = NMRD(
        u=u,
        atom_group=H_H2O,
        type_analysis="inter_molecular",
        number_i=20)
    res_h2o_inter = nmr_h2o_inter.run_analysis()

    nmr_peg_inter = NMRD(
        u=u,
        atom_group=H_PEG,
        type_analysis="inter_molecular",
        number_i=20)
    res_peg_inter = nmr_peg_inter.run_analysis()

The intermolecular contribution is
computed only between molecules belonging to the same chemical species.
For example, ``nmr_h2o_inter`` includes interactions between different
water molecules, but not between water and PEG molecules.

The water-PEG intermolecular contribution can be computed by selecting
water hydrogen atoms as the reference group (``atom_group``) and PEG
hydrogen atoms as the interacting partner group (``neighbor_group``):

.. code-block:: python

    nmr_h2o_peg = NMRD(
        u=u,
        atom_group=H_H2O,
        neighbor_group=H_PEG,
        number_i=20)
    res_h2o_peg = nmr_h2o_peg.run_analysis()

In this case, the analysis is already restricted to intermolecular
interactions between two different molecular species. Therefore, it is
not necessary to explicitly set ``type_analysis="inter_molecular"``.

Comparing the calculated spectra reveals that the
intramolecular contribution is larger than the intermolecular one for
this system (:ref:`Fig. 2 <fig:nmr-tutorial2>`). More importantly, the two contributions exhibit distinct
frequency dependences because they originate from different molecular
motions. Intramolecular relaxation is mainly governed by rotational motion and
internal molecular flexibility, whereas intermolecular relaxation reflects
translational diffusion and inter-molecular collisions.

.. _fig:nmr-tutorial2:

.. container:: figure

    .. image:: tutorial/nmr-intra-dm.png
        :class: only-dark
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. image:: tutorial/nmr-intra.png
        :class: only-light
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. container:: figurelegend

        Figure 2: Intramolecular :math:`^1\text{H}`-NMR
        relaxation rates :math:`R_{1 \text{R}}` (A) and
        Intermolecular :math:`^1\text{H}`-NMR relaxation
        rates :math:`R_{1 \text{T}}` (B) as a
        function of the frequency :math:`f` for 
        PEG and :math:`\text{H}_2\text{O}` separately.
        Results are shown for :math:`n_i = 0`.
