
Theory
======

Here, there theory behind the ``NMRDfromMD`` package is presented.

The system of interest is an ensemble of *identical* spins characterized
by a gyromagnetic ratio :math:`\gamma_I` and spin quantum number
:math:`I`. For :math:`^{1} \text{H}`, the most abundant isotope of hydrogen,
:math:`I = 1/2` and :math:`\gamma_I = 26.752` rad/T/s. For :math:`^{13} \text{C}`,
a natural and stable isotope of carbon, :math:`I = 1/2` and
:math:`\gamma_I = 6.728` rad/T/s :cite:`kowalewskiNuclearSpinRelaxation2006`.

.. 
    S.G.: Find a better place for that, there could be a paragraph explaining
    the limitation and hypothesis behind the thoery:
    - dipolar vs quadrupolar --> limited to dipolar because...
    - limited to systems where cross-correlation can be neglected --> what error?
    what systems?

One assumption behind the theory presented here is that cross-correlation terms
can be neglected; see Ref. :cite:`lippensT1RelaxationTime1993`.

When spin-lattice relaxation is dominated by fluctuations of the magnetic
dipole-dipole interaction, as is the case for protons in molecular systems,
the rates :math:`R_1 (\omega)` and :math:`R_2 (\omega)` are related to the
spectral densities :math:`J^{(m)}(\omega)` of these fluctuations via the
Bloembergen-Purcell-Pound (BPP) equations
:cite:`bloembergenRelaxationEffectsNuclear1948`:

.. math::
    :label: eq_BPP

    R_1 (\omega) & = & K \left[ J^{(1)} (\omega) + J^{(2)} (2 \omega) \right],

    R_2 (\omega) & = & K \left[ J^{(0)} (0) + 10 J^{(1)} (\omega)
    + J^{(2)} (2 \omega) \right] / 4,

where

.. math::

    K = \dfrac{3}{2}\left(\dfrac{\mu_0}{4 \pi}\right)^2 \hbar^2 \gamma^4 I (I+1),

where :math:`\mu_0` is the vacuum permeability, and :math:`\hbar` the Planck
constant (divided by :math:`2 \pi`). The constant :math:`K` has units of
:math:`\text{m}^6/\text{s}^2`. The spectral densities :math:`J^{(m)} (\omega)` in
Eq. :eq:`eq_BPP` can be obtained as the Fourier transform of the
autocorrelation functions :math:`G^{(m)}(\tau)`:

.. math::

    J^{(m)} (\omega) = \int_0^\infty G^{(m)} (\tau) \cos(\omega \tau) \mathrm d \tau.

The spectral densities are a measure of the distribution of the fluctuations
of :math:`G^{(m)}(\tau)` among different frequencies. They provide information
on the distribution of the power available for causing spin transitions at
different frequencies. The autocorrelation functions :math:`G^{(m)}(\tau)`
are given by

.. math::

    G^{(m)} (\tau) = \left< F_2^{(m)} [\textbf{r}_{ij} (t)]
    F_2^{*(m)} [\textbf{r}_{ij} (0)] \right>

where :math:`F_2^{(m)}` are complex functions of the vector
:math:`\textbf{r}_{ij}` between spin pairs, with norm :math:`r_{ij}` and
orientation :math:`\Omega_{ij}` with respect to a reference applied magnetic
field, assumed to be in the :math:`\textbf{e}_z` direction. The functions
:math:`F_2^{(m)}` are defined as

.. math::

    F_2^{(m)} [\textbf{r}_{ij} (t)] = \alpha_m \dfrac{Y_2^{(m)} [\Omega_{ij} (t)]}{r_{ij}^3 (t)}

where :math:`Y_2^{(m)}` are normalized spherical harmonics and
:math:`\alpha_0^2 = 16 \pi /5`, :math:`\alpha_1^2 = 8 \pi /15`, and
:math:`\alpha_2^2 = 32 \pi / 15`. Therefore, the correlation functions can be
written as:

.. math::

    G^{(m)} (\tau) = \dfrac{\alpha_m^2}{N} \sum_i \sum_{j \ne i}
    \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)} 
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)},

where :math:`N` is the number of spins.

Intra/inter contributions
-------------------------

Intra-molecular and inter-molecular contributions to :math:`R_1` and
:math:`R_2` can be extracted separately by splitting the correlation functions
as:

.. math::
    :label: G_intra

    G^{(m)}_\text{intra} (t) = \dfrac{\alpha_m^2}{N}
    \sum_i \sum_{j \in M_i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)}
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)},

.. math::
    :label: G_inter

    G^{(m)}_\text{inter} (t) = \dfrac{\alpha_m^2}{N}
    \sum_i \sum_{j \notin M_i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)}
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)},


where :math:`j \in M_i` and :math:`j \notin M_i` refer to spins from the same
molecule as :math:`i`, and from different molecules than :math:`i`,
respectively.

Intra-molecular relaxation is usually attributed to the rotational motion of
the molecules, and inter-molecular relaxation to their translational motion.
Although this assumption facilitates interpretation, it is not exact and
must be applied cautiously :cite:`hubbardTheoryNuclearMagnetic1963`.

Isotropic system
----------------

For isotropic systems, the correlation functions are proportional to each
other: :math:`G^{(0)} = 6 G^{(1)}`, and :math:`G^{(0)} = 6 / 4 G^{(2)}`
:cite:`becherMolecularDynamicsSimulations2021`. Therefore, there is no need to
calculate all three correlation functions, and :math:`G^{(0)} (t)` is usually
the only one computed, which considerably reduces the computational effort.

In that case, the rates :math:`R_1 (\omega)` and :math:`R_2 (\omega)` can be
written as:

.. math::

    R_1 &=&  \frac{K}{6} \left[ J^{(0)} (\omega_0) + 4 J^{(0)} (2 \omega_0) \right],

    R_2 &=& \frac{K}{6} \left[ J^{(0)} (0) + \frac{5}{2} J^{(0)} (\omega_0) + J^{(0)} (2 \omega_0) \right],

where

.. math::
    :label: F_2_0

    F_2^{(0)} [\textbf{r}_{ij} (t)] & = & \alpha_m \dfrac{Y_2^{(0)} [\Omega_{ij} (t)]}{r_{ij}^3 (t)}

    & = & \dfrac{3 \cos^2 \theta_\text{ij} (t) - 1}{r_{ij}^3 (t)}

Here, we check the validity of the relation
:math:`G^{(0)} = 6 G^{(1)} = 6 / 4 G^{(2)}` on a simple bulk water system with
4000 molecules, similar to the approach taken in
:cite:`becherMolecularDynamicsSimulations2021` with glycerol. The proportionality
relation is well verified (Figure below).

.. image:: ../figures/illustrations/bulk-water/effect_of_anisotropy-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/illustrations/bulk-water/effect_of_anisotropy-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: Test of the validity of the relation
    :math:`G^{(0)} = 6 G^{(1)} = 6 / 4 G^{(2)}` on a bulk water system; see text
    for details.