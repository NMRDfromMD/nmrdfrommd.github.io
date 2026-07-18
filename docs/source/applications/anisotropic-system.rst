.. include:: ../additional/links.rst
.. _anisotropic-label:

Nanoconfined water
==================

Here, the :math:`^1\text{H}` NMR relaxation rate :math:`R_1` of water confined
within a silica nanoslit is calculated from molecular dynamics simulations. The
temperature is :math:`T = 300~\mathrm{K}`, and the water model is
:math:`\text{TIP4P}-\epsilon`. Simulation details are provided in
:ref:`simulation-methods`. For such an anisotropic system, all three
correlation functions, :math:`G^{(1)}`, :math:`G^{(2)}`, and
:math:`G^{(3)}`, must be evaluated. Hydrogen atoms belonging to both the water
molecules and the hydroxyl groups on the silica surface contribute to the
relaxation, particularly because water forms hydrogen bonds with the surface.

Using NMRDfromMD, we calculated both the total relaxation rate and the
contribution arising from water-surface interactions. Our results indicate
that the contribution from water-silica interactions is approximately one
order of magnitude smaller than the total relaxation rate, indicating that
most of the relaxation originates from the water itself through its rotational
motion and intermolecular interactions (:ref:`Fig. 1 <fig:confined>`). This is
mainly due to the relatively small number of hydrogen atoms on the silica
surface (92), compared with the 1204 hydrogen atoms belonging to the water
molecules.

.. _fig:confined:

.. container:: figure

    .. image:: anisotropic-system/nmr-water-silica-dm.png
        :class: only-dark
        :alt: NMR results obtained from the GROMACS simulation of water in silica

    .. image:: anisotropic-system/nmr-water-silica.png
        :class: only-light
        :alt: NMR results obtained from the GROMACS simulation of water in silica

    .. container:: figurelegend

        Figure 1. (A) :math:`^1\text{H}` NMR relaxation rate :math:`R_1` of
        water confined within a silica slit. The contribution arising from
        water-silica interactions is shown by the pink pentagons. The dashed
        line marks the expected relaxation rate of
        :math:`0.28~\mathrm{s}^{-1}` for the same TIP4P water model at the
        same temperature, :math:`T = 300~\mathrm{K}`. (B) Snapshot of the
        molecular dynamics system, with water molecules shown in red and
        white, silicon atoms in yellow, and sodium counterions in blue.

Although the direct contribution of water-surface interactions to NMR
relaxation is small, confinement still significantly affects the relaxation
rate by modifying the structure and dynamics of the confined water. In
particular, confinement generally slows the molecular dynamics of water,
resulting in altered relaxation rates. Furthermore, the confined water can be
viewed as consisting of two populations: a slow, highly structured
interfacial population located near the silica surface, and a bulk-like
population farther from the surface whose properties are only weakly affected
by confinement :cite:`gravelleIntermittentMolecularMotion2025`. As a
consequence, the total relaxation rate is slightly larger than that measured
for bulk water under the same conditions (:ref:`Fig. 1 <fig:confined>`).
