.. _lysozyme-label:

Lysozyme in water
=================

The system consists of a lysozyme (HEWL) and 594 water molecules. The
temperature is :math:`T = 300~\mathrm{K}`, and the water model is
:math:`\text{TIP4P}-\epsilon`. Simulation details are provided in
:ref:`simulation-methods`.

The total NMR relaxation rate :math:`R_1` exhibits a strong frequency
dependence across the entire accessible frequency range, in contrast to simple
bulk liquids, where :math:`R_1` reaches a low-frequency plateau. Decomposing
the signal into water and lysozyme contributions reveals that the
low-frequency dispersion is dominated by the protein
(:ref:`Figure 1 <fig:lysozyme>`). This reflects the slow rotational tumbling of
the lysozyme molecule, which has a much longer correlation time than the small
solvent molecules.

The water contribution alone also shows a residual frequency dependence at low
frequencies, which is absent in bulk water. In bulk water,
:math:`R_1` reaches a plateau below approximately
:math:`2 \cdot 10^3\,\text{MHz}`, indicating that molecular motion is fast
relative to the NMR timescale. In contrast, water molecules in contact with
the lysozyme surface show dispersion extending down to approximately
:math:`10\,\text{MHz}`. This is consistent with the expected slowdown of the
translational and rotational dynamics of adsorbed water molecules, whose
motion is partially constrained by interactions with the protein surface.

.. _fig:lysozyme:

.. container:: figure
        
    .. image:: lysozyme-in-water/nmr-water-hewl-dm.png
        :class: only-dark
        :alt: NMR results obtained from the LAMMPS simulation of water and lysozyme

    .. image:: lysozyme-in-water/nmr-water-hewl.png
        :class: only-light
        :alt: NMR results obtained from the LAMMPS simulation of water and lysozyme

    .. container:: figurelegend

        Figure 1. (A) NMR relaxation rate :math:`R_1` for the
        lysozyme--water system. The spectra for water alone and lysozyme alone
        are shown, alongside the contribution for the water-lysozyme interaction.
        (B) Snapshot of the molecular dynamics simulation.
