.. _correlation-spectra:


Lorentzian and diffusion models
===============================

This section describes how molecular dynamics trajectories are used to
identify the appropriate models for NMR relaxation. The dipolar
fluctuations extracted from the trajectories are first converted into
correlation functions and spectral densities. Their functional forms
then reveal the underlying relaxation mechanisms: a Lorentzian model
associated with rotational diffusion and a diffusion-based model
associated with translational motion.

Although the behaviour of a single pair of nuclei is informative, NMR
relaxation is a collective property. The relaxation rates are obtained
by averaging the fluctuations of :math:`F_2^{(0)}` over all relevant
pairs of nuclei and over the complete molecular dynamics trajectory.

From the fluctuating quantities :math:`F_2^{(0)}` summed over all
available pairs of spins, one can extract the two correlation functions
:math:`G_\textrm{intra}^{(0)}` and :math:`G_\textrm{inter}^{(0)}`
(see Eqs. :eq:`G_intra` and :eq:`G_inter`).

For comparison, the results obtained at two different temperatures,
275 and 300 K, are reported.

At short times, :math:`t < 40` ps, the intramolecular correlation
functions follow an exponential decay,

.. math::
    :label: eq_exp_G

    G_\text{intra} (t) =
    G_\text{intra} (0) \exp{(-t / \tau_\text{intra})},

where :math:`\tau_\text{intra} = 6.3` ps was used for :math:`T = 300`
K and :math:`\tau_\text{intra} = 3.2` ps was used for :math:`T = 275`
K, see the figure below.

An exponential decay, such as Eq. :eq:`eq_exp_G`, is characteristic of
processes governed by a single correlation time. Such behaviour is
commonly used to describe isotropic rotational diffusion and provides a
good approximation for the short-time intramolecular dynamics of liquid
water :cite:`lippensT1RelaxationTime1993`.

The intermolecular correlation functions, however, display a different
behaviour. They follow an exponential decay only at short times (a few
tens of picoseconds). At longer times, translational diffusion
continually brings new molecular neighbours into and out of the local
environment. This leads to the characteristic power-law decay
:math:`G_\text{inter}(t) \sim t^{-3/2}`,
which is a hallmark of diffusive dynamics.

This scaling has been predicted theoretically for freely diffusing
particles, and analytical expressions were derived by Ayant et al.
:cite:`ayantCalculDensitesSpectrales1975` and Hwang and Freed
:cite:`hwangDynamicEffectsPair2008` in the context of hard-sphere
diffusion.

Following Ref. :cite:`grivetNMRRelaxationParameters2005`, we refer to
this description as ADHF.

While the ADHF description captures the long-time decay of the
correlation function, it does not lead to a Lorentzian spectral
density. This reflects the fact that intermolecular relaxation is not
controlled by a single correlation time, but by translational diffusion
and molecular exchange processes.

.. image:: correlation-spectra/spectra-dm.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: correlation-spectra/spectra.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: (A) Intramolecular dipolar correlation function
    :math:`G^{(0)}_\mathrm{intra}(t)`. (B) Intermolecular dipolar
    correlation function :math:`G^{(0)}_\mathrm{inter}(t)`, shown on
    a log-log scale. The dashed line indicates the long-time scaling
    :math:`G(t) \propto t^{-3/2}`, characteristic of translational
    diffusion. (C) Intramolecular spectral density
    :math:`J^{(0)}_\mathrm{intra}(f)`. The solid line corresponds to a
    Lorentzian fit based on a single correlation time approximation.
    (D) Intermolecular spectral density
    :math:`J^{(0)}_\mathrm{inter}(f)`. The solid line shows the
    analytical prediction based on adsorption-diffusion.

The different physical origins of intra- and intermolecular relaxation
are reflected in their spectral densities. Intramolecular relaxation is
dominated by molecular rotation and can therefore be described by a
single correlation time approximation. In contrast, intermolecular
relaxation is affected by translational diffusion and requires a
different description.

The intramolecular spectrum :math:`J_\textrm{intra}^{(0)}` can be
reasonably well described by a Lorentzian:

.. math::
    :label: eq_lorenzian_G

    J_\text{intra} (f) =
    G_\text{intra} (0)
    \dfrac{2 \tau_\text{c}}
    {1 + 4 \pi^2 f^2 \tau_\text{c}^2}

using :math:`\tau_\text{c} = 6.3` ps and
:math:`G(0) = 56300~\mathrm{Å}^{-6}\,\mathrm{ps}^{-2}` for
:math:`T = 300` K, and
:math:`\tau_\text{c} = 3.2` ps and
:math:`G(0) = 59500~\mathrm{Å}^{-6}\,\mathrm{ps}^{-2}` for
:math:`T = 275` K.

The intermolecular spectral density
:math:`J_\mathrm{inter}^{(0)}`, however, does not exhibit a Lorentzian
plateau at low frequencies. This deviation is directly related to the
long-time power-law decay observed in the corresponding correlation
function :math:`G_\mathrm{inter}^{(0)}(t)`.

In particular, the algebraic decay
:math:`G(t) \sim t^{-3/2}` is characteristic of translational diffusion
and indicates that the relaxation dynamics is not governed by a single
characteristic correlation time.

In this regime, and following closely
Ref. :cite:`gravelleAdsorptionKineticsOpen2019`, the spectral density
can be described by adsorption-diffusion, which yields an analytical
expression based on first-passage statistics in a diffusive reservoir:

.. math::
    :label: eq_spectrum_sqrt

    J_\mathrm{inter}(f)
    = \left[1 + A + B \sqrt{2 \pi f}\right]^{-1}.

Within this framework, the parameters have a direct physical
interpretation:
:math:`A = k r / D` and :math:`B = r / \sqrt{D}`,
where :math:`r` is the molecular radius, :math:`D` the self-diffusion
coefficient, and :math:`k` a phenomenological exchange rate
(with units of m/s) describing effective adsorption-desorption kinetics.

The resulting frequency dependence captures the crossover from a
quasi-plateau at low frequency to a diffusion-dominated decay at higher
frequencies.

As shown in panel (D), the ADHF model provides a good description of
the molecular dynamics results up to approximately
:math:`10^4\,\mathrm{MHz}`.
