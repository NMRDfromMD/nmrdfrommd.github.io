.. _bulk-water-label:

Microscopic mechanisms of relaxation
====================================

This section illustrates how the microscopic motions observed in a
molecular dynamics trajectory give rise to the macroscopic
:math:`^1\mathrm{H}`-NMR relaxation rates measured experimentally.
Using bulk liquid water as a simple example, we follow the complete
chain of quantities entering the relaxation calculation, from the
motion of individual pairs of nuclei to the dipolar correlation
functions and, finally, to the relaxation spectra.

Throughout this section, intramolecular and intermolecular
contributions are analysed separately, highlighting how rotational and
translational molecular motions contribute differently to NMR
relaxation.

MD system
---------

.. image:: ../figures/tutorials/bulk-water/water-dark-square.png
    :class: only-dark
    :alt: Water molecules simulated with lammps - NMR relaxation time calculation
    :width: 250
    :align: right

.. image:: ../figures/tutorials/bulk-water/water-light-square.png
    :class: only-light
    :alt: Water molecules simulated with lammps - NMR relaxation time calculation
    :width: 250
    :align: right

The system consists of a bulk liquid containing :math:`N = 2000` water
molecules. The simulation box was cubic, with equilibrium dimensions of
:math:`(X\,\text{nm})^3`. Simulations were performed using the open-source
code LAMMPS. The system was first equilibrated in the NPT ensemble at a
temperature of :math:`T = 300\,\text{K}` and a pressure of
:math:`p = 1\,\text{atm}`. The trajectory was then recorded during a
:math:`10\,\text{ns}` production run performed in the NVT ensemble using a
timestep of :math:`2\,\text{fs}`. The atomic positions were written to the
:file:`prod.xtc` trajectory file every :math:`\Delta t = 2\,\text{ps}`.

In the following, we progressively connect molecular motion in the
trajectory to correlation functions and spectral densities, in order
to identify the distinct microscopic mechanisms contributing to
NMR relaxation.

Isotropic consistency of :math:`G^{(m)}` functions
--------------------------------------------------

Before analysing the molecular dynamics, we first verify one of the
central assumptions used throughout the package. For isotropic bulk
liquids, the three second-rank dipolar correlation functions are
expected to satisfy :cite:`hubbardTheoryNuclearMagnetic1963`

.. math::
    
    G^{(0)} = 6 G^{(1)} = \frac{6}{4} G^{(2)}.

For an isotropic bulk liquid, no direction in space is preferred.
The three correlation functions :math:`G^{(m)}` differ only by their
spherical harmonic order :math:`m = 0, 1, 2`. Because all orientations
are equally probable, the orientation average cannot depend on :math:`m`,
and the functions :math:`G^{(m)}` are therefore proportional to one
another. The numerical prefactors :math:`1, 6, 6/4` arise solely from
the explicit forms of the rank-2 spherical harmonics :math:`Y_2^m`,
whose squared moduli satisfy
:math:`|Y_2^0|^2 : |Y_2^1|^2 : |Y_2^2|^2 = 1 : 3 : 3`, combined with
the symmetry factor accounting for :math:`m` > 0 pairs.

Verifying this proportionality provides a simple consistency check
before analysing the relaxation spectra. A similar bencmark was done, for instance
in Ref. :cite:`becherMolecularDynamicsSimulations2021` with glycerol. Our results
show that the proportionality relation is well verified (See figure below).

.. image:: microscopic-origin/proportionality-G-dm.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: microscopic-origin/proportionality-G.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: A) Comparison between :math:`G^{(0)}`, :math:`6 G^{(1)}` and :math:`\frac{6}{4} G^{(2)}`
    obtained from a bulk water molecular dynamics. B) Relaxation rate, :math:`R_1`,
    inter-molecular relaxation rate, :math:`R_\text{1 T}`, and intra-molecular relaxation rate, :math:`R_\text{1 R}`, 
    as a function of the frequency :math:`f`. The temperature is 300 K, and 
    the total number of water molecules is 3000.

Our results also show that the total relaxation is dominated by intra-molecular contribution,
as expected for water under ambient conditions :cite:`singerMolecularDynamicsSimulations2017`.

To understand the origin of these relaxation spectra, it is useful to
return to the atomic trajectories themselves. The dipolar interaction
between two nuclear spins depends on both their separation
:math:`r_{ij}` and their relative orientation
:math:`\Omega_{ij}`. Following the evolution of these quantities in
time therefore provides direct insight into the microscopic origin of
NMR relaxation.

Given that the system is isotropic, the correlation functions
satisfy :math:`G^{(0)} = 6 G^{(1)} = 6/4 G^{(2)}`. It is therefore
sufficient to compute only :math:`G^{(0)}` and the corresponding
spectral density :math:`J^{(0)}`. In this case, the angular dependence
reduces to the polar angle :math:`\theta_{ij}`, since
:math:`Y_2^{(0)}` is independent of the azimuthal angle
:math:`\varphi`.

Intra-molecular contribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first consider the two hydrogen atoms belonging to the same water
molecule. Because the present water model (TIP4P/:math:`\epsilon`) is rigid,
the internuclear distance :math:`r_{ij}` remains essentially constant
(within the tolerance of the SHAKE algorithm used to enforce molecular
rigidity). The only quantity that changes with time is therefore the
orientation of the H-H vector with respect to the external magnetic field,
described by the polar angle :math:`\theta_{ij}` (see the figure below).

The dipolar interaction entering the relaxation equations is described by
:math:`F_2^{(0)}` (Eq. :eq:`F_2_0`), which depends on both the internuclear
distance and its orientation. Since :math:`r_{ij}` is nearly constant for a
rigid water molecule, the fluctuations of :math:`F_2^{(0)}` arise almost
entirely from the rotational motion of the molecule through changes in
:math:`\theta_{ij}`.

For a typical H-H distance :math:`a \approx 1.51\,\mathrm{Å}`,
:math:`F_2^{(0)}` varies between

.. math::

    \frac{3\cos^2 0 - 1}{a^3}
    \approx 0.58~\mathrm{Å}^{-3}

and

.. math::

    \frac{3\cos^2 (\pi/2) - 1}{a^3}
    \approx -0.29~\mathrm{Å}^{-3},

corresponding to the H--H vector being parallel and perpendicular to the
magnetic field, respectively. Thus, although the internuclear distance is
essentially fixed, molecular rotation produces substantial fluctuations of
the dipolar interaction.

.. image:: microscopic-origin/intra-molecular-dm.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: microscopic-origin/intra-molecular.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: A) :math:`\theta_{ij}` as a function of the time :math:`t`, where :math:`i` and :math:`j`
    refer to two hydrogen atoms located within the same water molecule. B) :math:`r_{ij}` as a function of 
    time. C) :math:`F_{2}^{(0)}` as a function of time. The temperature is 300 K, and 
    the total number of water molecules is 4000.

Inter-molecular contribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We now consider two hydrogen atoms belonging to different water
molecules. In contrast with the intramolecular case, both the
internuclear distance and the relative orientation fluctuate because
of translational diffusion. In that case, :math:`r_{ij}` fluctuates significantly between :math:`\approx 2.5 A`,
corresponding to molecules occupying their respective hydration shell, to larger values
(potentially as large as the box permits). As can be seen, the function :math:`F_{0}^{(2)}`
reaches its largest absolute values, here about :math:`0.02 \mathrm{Å}^{-3}`,
when :math:`r_{ij}` is the shorter.

.. image:: microscopic-origin/inter-molecular-dm.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: microscopic-origin/inter-molecular.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: A) :math:`\theta_{ij}` as a function of the time :math:`t`, where :math:`i` and :math:`j`
    refer to two hydrogen atoms located within two different water molecules. B) :math:`r_{ij}` as a function of 
    time. C) :math:`F_{2}^{(0)}` as a function of time. The temperature is 300 K, and 
    the total number of water molecules is 4000.

From pair dynamics to collective relaxation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although the behaviour of a single pair of nuclei is informative, NMR
relaxation is a collective property. The relaxation rates are obtained
by averaging the fluctuations of :math:`F_2^{(0)}` over all relevant
pairs of nuclei and over the complete molecular dynamics trajectory.

From the fluctuating quantities :math:`F_{0}^{(2)}` summed up over all the available pair of 
spins, one can extract the two correlation functions :math:`G_\textrm{intra}^{(0)}` and
:math:`G_\textrm{inter}^{(0)}` (see Eqs. :eq:`G_intra` and :eq:`G_inter`). For comparison,
the results obtained with two different temperatures 275 and 300 K are reported.

At short time :math:`t < 40` ps, the intra-molecular correlation functions follow and
a decreasing exponential,

.. math::
    :label: eq_exp_G

    G_\text{intra} (t) = G_\text{intra} (0)  \exp{(-t / \tau_\text{intra})},

where :math:`\tau_\text{intra} = 6.3` ps was used for :math:`T = 300` K 
and :math:`\tau_\text{intra} = 3.2` ps was used for :math:`T = 275` K, see the figure 
below. An exponential decay, such as Eq. :eq:`eq_exp_G`, is common for process governed by a single
characteristic correlation time. Such behaviour is commonly
used to describe isotropic rotational diffusion and provides a good
approximation for the short-time intramolecular dynamics of liquid
water :cite:`lippensT1RelaxationTime1993`.

The inter-molecular correlation functions, however, scale as an
exponential [i.e. Eq. :eq:`eq_exp_G`] only for time shorter than a 
few tens of pico-second, and show a clear scaling as :math:`G_\text{inter} (t) \sim t^{-3/2}`
for large time which is a characteristic signature of the diffusion
process controlling the motion of the molecules. 

At long times, translational diffusion continually brings new
molecular neighbours into and out of the local environment. This
leads to a characteristic long-time decay
:math:`G_\text{inter}(t) \sim t^{-3/2}`, which is a hallmark of
diffusive dynamics. This scaling has been predicted theoretically for freely diffusing
particles, and analytical expressions were derived by Ayant et al. :cite:`ayantCalculDensitesSpectrales1975`
and Hwang and Freed :cite:`hwangDynamicEffectsPair2008` in the context of hard-sphere diffusion.
Following Ref. :cite:`grivetNMRRelaxationParameters2005`, we refer
to this description as ADHF.

While the ADHF description captures the long-time decay of the
correlation function, it does not lead to a Lorentzian spectral
density. This motivates the use of alternative theoretical
descriptions for the low-frequency behaviour of
:math:`J_\text{inter}(f)`

.. image:: microscopic-origin/spectra-dm.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: microscopic-origin/spectra.png
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: (A) Intramolecular dipolar correlation function :math:`G^{(0)}_\mathrm{intra}(t)`.
    (B) Intermolecular dipolar correlation function :math:`G^{(0)}_\mathrm{inter}(t)`,
    shown on a log-log scale. The dashed line indicates the long-time scaling
    :math:`G(t) \propto t^{-3/2}`, characteristic of translational diffusion.
    (C) Intramolecular spectral density :math:`J^{(0)}_\mathrm{intra}(f)`.
    The solid line corresponds to a Lorentzian fit based on a single correlation time
    approximation.
    (D) Intermolecular spectral density :math:`J^{(0)}_\mathrm{inter}(f)`.
    The solid line shows the analytical ADHF prediction capturing diffusion-driven
    deviations from Lorentzian behavior at low frequencies.

The intra molecular spectrum :math:`J_\textrm{intra}^{(0)}` can be reasonably
well adjusted by a Lorentzian

.. math::
    :label: eq_lorenzian_G

    J_\text{intra} (f) = G_\text{intra} (0) \dfrac{2 \tau_\text{c}}{1 + 4 \pi^2 f^2 \tau_\text{c}^2}

using :math:`\tau_\text{c} = 6.3` ps and :math:`G(0) = 56300` A⁻⁶ ps⁻² for :math:`T = 300` K
and :math:`\tau_\text{c} = 3.2` ps and :math:`G(0) = 59500` A⁻⁶ ps⁻² for :math:`T = 275` K. 

The intermolecular spectral density :math:`J_\mathrm{inter}^{(0)}`, however, does not follow a Lorentzian
plateau at low frequencies. This deviation is directly related to the long-time
power-law decay observed in the corresponding correlation function
:math:`G_\mathrm{inter}^{(0)}(t)`.

In particular, the algebraic decay :math:`G(t) \sim t^{-3/2}` is caracteristic of
translational diffusion and indicates that the relaxation dynamics is not governed
by a single characteristic correlation time. In this regime, and following closely
Ref. :cite:`gravelleAdsorptionKineticsOpen2019`, the spectral density can be
described by the adsorption-diffusion hard-sphere (ADHF) model, which yields an
analytical expression based on first-passage statistics in a diffusive reservoir:

.. math::
    :label: eq_spectrum_sqrt

    J_\mathrm{ADHF}(f) = \left[ 1 + A + B \sqrt{2 \pi f} \right]^{-1}.

Within this framework, the parameters have a direct physical interpretation:
:math:`A = k r / D` and :math:`B = r / \sqrt{D}`, where :math:`r` is the molecular
radius, :math:`D` the self-diffusion coefficient, and :math:`k` a phenomenological
exchange rate (m·s⁻¹) describing effective adsorption-desorption kinetics.

The resulting frequency dependence captures the crossover from a quasi-plateau at
low frequency to a diffusion-dominated decay at higher frequencies. As shown in
panel (D), the ADHF model provides a good description of the molecular dynamics
results up to approximately :math:`10^4\,\mathrm{MHz}`.
