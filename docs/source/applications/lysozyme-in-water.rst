.. _lysozyme-label:

Lysozyme in water
=================

MD system
---------

.. image:: ../figures/illustrations/lysozyme-in-water/snapshot-dark.png
    :class: only-dark
    :alt: lysozyme in water simulated with GROMACS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

.. image:: ../figures/illustrations/lysozyme-in-water/snapshot-light.png
    :class: only-light
    :alt: lysozyme in water simulated with GROMACS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

The system is made of a lysozyme (HEWL) with 594 water molecules, which
corresponds to water-to-protein mass ratio of :math:`73\,\%`.
The simulation was made using GROMACS using a timestep of :math:`2\,\text{fs}`.
The simulation was performed using GROMACS with a :math:`2\,\text{fs}`. timestep,
a :math:`100\,\text{ns}`, and trajectories recorded every 1 ps at 300 K.

.. admonition:: Note
    :class: non-title-info

    If you are new to NMRDfromMD, it is recommended to begin with the :ref:`isotropic-label`.

Results
-------

The total NMR relaxation rate :math:`R_1` exhibits a strong frequency dependence
across the entire accessible frequency range, in contrast to simple bulk liquids
where :math:`R_1` reaches a low-frequency plateau. Decomposing the signal into
water and lysozyme contributions reveals that the low-frequency dispersion is
dominated by the protein. This reflects the slow rotational tumbling of the
lysozyme molecule, which has a much longer correlation time than small solvent
molecules.

The water contribution alone also shows a residual frequency dependence at low
frequencies, which is absent in pure bulk water. In bulk water, :math:`R_1`
reaches a plateau below approximately :math:`2 \cdot 10^3\,\text{MHz}`,
indicating that molecular motion is fast relative to the NMR timescale. In
contrast, water molecules in contact with the lysozyme surface show dispersion
extending down to approximately :math:`10\,\text{MHz}`. This is consistent with
the expected slowdown of translational and rotational dynamics of adsorbed water
molecules, whose motion is partially constrained by interactions with the protein
surface.

.. image:: ../figures/illustrations/lysozyme-in-water/R1_spectra-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water and lysozyme

.. image:: ../figures/illustrations/lysozyme-in-water/R1_spectra-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water and lysozyme

.. container:: figurelegend

    Figure: NMR relaxation rate :math:`R_1` for the lysozyme-water system.
    The spectra for water alone and lysozyme alone are also given.

This comparison illustrates a key advantage of MD-based NMR relaxation
calculations: the ability to decompose the total signal into contributions from
distinct molecular species, providing physical insight that is inaccessible from
experiment alone.

.. image:: ../figures/illustrations/lysozyme-in-water/R1_spectra_water-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water and lysozyme

.. image:: ../figures/illustrations/lysozyme-in-water/R1_spectra_water-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water and lysozyme

.. container:: figurelegend

    Figure: Contribution to the NMR relaxation rate :math:`R_1` from the water only,
    comparing the residual water in contact with the lysozyme, and pure bulk water.
