Mathematical framework
======================

The system of interest here is an ensemble of identical spins characterized
by a gyromagnetic ratio :math:`\gamma_I` and spin quantum number
:math:`I`. For :math:`^{1} \text{H}`, the most abundant isotope of hydrogen,
:math:`I = 1/2` and :math:`\gamma_I = 26.752` :math:`\text{rad} \, \text{T}^{-1} \, \text{s}^{-1}`. For :math:`^{13} \text{C}`,
a natural and stable isotope of carbon, :math:`I = 1/2` and
:math:`\gamma_I = 6.728` :math:`\text{rad} ~ \text{T}^{-1} \, \text{s}^{-1}` :cite:`kowalewskiNuclearSpinRelaxation2006`.

One assumption behind the theory presented here is that cross-correlation terms
can be neglected, i.e. that correlated fluctuations of the dipolar interactions
do not contribute significantly to the relaxation. This approximation is
discussed in Ref. :cite:`lippensT1RelaxationTime1993`.

When spin-lattice relaxation is dominated by fluctuations of the magnetic dipole-dipole
interaction, which is often the dominant mechanism for protons in molecular systems,
the relaxation rates :math:`R_1(\omega_0)` and :math:`R_2(\omega_0)` can be expressed
in terms of the spectral densities :math:`J^{(m)}(\omega)` through the Bloembergen-Purcell-Pound (BPP)
formalism :cite:`bloembergenRelaxationEffectsNuclear1948`. Here, :math:`\omega_0 = \gamma_I B_0`
denotes the Larmor angular frequency in the static magnetic field :math:`B_0`.
:math:`R_1` describes longitudinal (spin-lattice) relaxation, governing recovery
of the magnetization along the static field, whereas :math:`R_2` describes transverse
(spin-spin) relaxation, governing the decay of phase coherence in the plane perpendicular to the field.
Within the BPP formalism, the longitudinal relaxation rate is given by

.. math::
    :label: eq_BPP_R1

    R_1(\omega_0) = K \left[ J^{(1)}(\omega_0) + J^{(2)}(2\omega_0) \right]

while the transverse relaxation rate is given by

.. math::
    :label: eq_BPP_R2

    R_2(\omega_0) = K \left[ J^{(0)}(0) + 10 J^{(1)}(\omega_0) + J^{(2)}(2\omega_0) \right] / 4

where

.. math::

    K = \dfrac{3}{2}\left(\dfrac{\mu_0}{4 \pi}\right)^2 \hbar^2 \gamma^4 I (I+1),

where :math:`\mu_0` is the vacuum permeability, and :math:`\hbar` the Planck
constant (divided by :math:`2 \pi`). The constant :math:`K` has units of
:math:`\text{m}^6/\text{s}^2`. The spectral densities :math:`J^{(m)} (\omega)` in
Eq. :eq:`eq_BPP` can be obtained as the Fourier transform of the
autocorrelation functions :math:`G^{(m)}(\tau)`:

.. math::

    J^{(m)} (\omega) = 2 \int_0^\infty G^{(m)} (\tau) \cos(\omega \tau) \mathrm d \tau.

The spectral densities are a measure of the distribution of the fluctuations
of :math:`G^{(m)}(\tau)` among different frequencies. They provide information
on the distribution of the power available for causing spin transitions at
different frequencies. The autocorrelation functions :math:`G^{(m)}(\tau)`
are given by

.. math::

    G^{(m)} (\tau) = \left< F_2^{(m)} [\textbf{r}_{ij} (\tau)]
    F_2^{*(m)} [\textbf{r}_{ij} (0)] \right>_{ij, \tau}

where :math:`F_2^{(m)}` are complex functions of the vector
:math:`\textbf{r}_{ij}` between spin pairs, with norm :math:`r_{ij}` and
orientation :math:`\Omega_{ij}` with respect to a reference applied magnetic
field, assumed to be in the :math:`\textbf{e}_z` direction. Here,
:math:`\left< \cdot \right>_{ij, \tau}` denotes an average over all spin
pairs and over time origins :math:`\tau`.

Note, some textbooks absorb numerical factors into :math:`J (m)`, leading to
different-looking BPP formulas [Eq. :eq:`eq_BPP`].

The functions :math:`F_2^{(m)}` are defined as

.. math::

    F_2^{(m)} [\textbf{r}_{ij} (t)] = \alpha_m \dfrac{Y_2^{(m)} [\Omega_{ij} (t)]}{r_{ij}^3 (t)}

where :math:`Y_2^{(m)}` are normalized spherical harmonics and
:math:`\alpha_0^2 = 16 \pi /5`, :math:`\alpha_1^2 = 8 \pi /15`, and
:math:`\alpha_2^2 = 32 \pi / 15`. Therefore, the correlation functions can be
written as:

.. math::

    G^{(m)} (\tau) = \dfrac{\alpha_m^2}{N} \sum_i \sum_{j \ne i}
    \left< \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)} 
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)} \right>_{\tau},

where :math:`N` is the number of spins, the double sum runs over all ordered
pairs :math:`(i,j)`, and the :math:`1/N` factor averages over spin origins
:math:`i`. Here, :math:`\left< \cdot \right>_{\tau}` denotes an average over all time origins.

Intra/inter contributions
-------------------------

Intra-molecular (R) and inter-molecular (T) contributions to :math:`R_1` and
:math:`R_2` can be extracted separately by splitting the correlation functions
as:

.. math::
    :label: G_intra

    G^{(m)}_\text{R} (\tau) = \dfrac{\alpha_m^2}{N}
    \left< \sum_i \sum_{j \in M_i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)}
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)} \right>_{\tau},

.. math::
    :label: G_inter

    G^{(m)}_\text{T} (\tau) = \dfrac{\alpha_m^2}{N}
    \left< \sum_i \sum_{j \notin M_i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)}
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)} \right>_{\tau},


where :math:`j \in M_i` and :math:`j \notin M_i` refer to spins from the same
molecule as :math:`i` (but different from :math:`i`), and from different molecules than :math:`i`,
respectively.

In isotropic liquids, the intra-molecular contribution is often dominated by
rotational motion, whereas the inter-molecular contribution is often dominated
by translational motion. This interpretation is approximate, however, as both
rotational and translational dynamics contribute to the time dependence of the
dipolar interaction, and the two cannot be fully decoupled
:cite:`hubbardTheoryNuclearMagnetic1963`.

Isotropic system
----------------

For isotropic systems, the correlation functions are proportional to each
other: :math:`G^{(0)} = 6 G^{(1)}`, and :math:`G^{(0)} = 6 / 4 G^{(2)}`
:cite:`becherMolecularDynamicsSimulations2021`. Therefore, there is no need to
calculate all three correlation functions, and :math:`G^{(0)} (\tau)` is usually
the only one computed, which considerably reduces the computational effort.

In that case, the rates :math:`R_1 (\omega_0)` and :math:`R_2 (\omega_0)` can be
written as:

.. math::
    :label: eq_BPP_R1_iso

    R_1 (\omega_0) =  K \left[ J^{(0)} (\omega_0) + 4 J^{(0)} (2 \omega_0) \right] / 6,

.. math::
    :label: eq_BPP_R2_iso

    R_2 (\omega_0) = K \left[ \frac{3}{2} J^{(0)} (0) + \frac{5}{2} J^{(0)} (\omega_0) + J^{(0)} (2 \omega_0) \right] / 6,

where

.. math::
    :label: F_2_0

    F_2^{(0)} [\textbf{r}_{ij} (t)] = \alpha_0 \dfrac{Y_2^{(0)} [\Omega_{ij} (t)]}{r_{ij}^3 (t)}
    = \dfrac{3 \cos^2 \theta_\text{ij} (t) - 1}{r_{ij}^3 (t)}.
