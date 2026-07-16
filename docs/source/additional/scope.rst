Scope
=====

Current capabilities
--------------------

NMRDfromMD is currently designed for the prediction and analysis of
:math:`^{1} \text{H}` nuclear magnetic relaxation dispersion (NMRD) profiles
from molecular dynamics trajectories. The current implementation focuses on
relaxation mechanisms dominated by magnetic dipole-dipole interactions, which
can be directly described through the time-dependent fluctuations of
internuclear vectors. In particular, :math:`^{1} \text{H}` relaxation can be
calculated from the rotational autocorrelation functions of relevant
interatomic vectors extracted from MD trajectories, followed by the
calculation of the corresponding spectral densities
:cite:`singerMolecularDynamicsSimulations2017, valiyaparambathuMolecularModesElucidate2024`.

The current framework is applicable to both isotropic and anisotropic
molecular systems. For isotropic systems, such as bulk liquids or freely
tumbling molecules, the relaxation analysis relies on the assumption of
isotropic rotational diffusion. For anisotropic systems, such as confined
liquids or molecules interacting with surfaces, the full
orientation-dependent dynamics of internuclear vectors are retained, allowing
the effect of restricted or heterogeneous molecular motions on the relaxation
dispersion profile to be investigated
:cite:`amaro-estradaQuantifyingEffectSpatial2022, gravelleIntermittentMolecularMotion2025`.
This enables the study of systems where deviations from simple isotropic rotational
diffusion play a key role in determining the observed NMR relaxation behaviour.

Future extensions
-----------------

Currently, NMRDfromMD is focused on the prediction of
:math:`^{1}\mathrm{H}` nuclear magnetic relaxation dispersion. This focus is
motivated by the relatively simple and well-established relationship between
proton relaxation rates and molecular motion, where relaxation is often
dominated by dipole-dipole interactions that can be directly described
through the autocorrelation functions of internuclear vectors extracted from
molecular dynamics trajectories.

Extending the framework to other nuclei, such as
:math:`^{13}\mathrm{C}`, :math:`^{15}\mathrm{N}`, and
:math:`^{19}\mathrm{F}`, is more challenging because the dominant relaxation
mechanisms depend strongly on the chemical environment. For example,
:math:`^{13}\mathrm{C}` relaxation in protonated carbons is often governed by
heteronuclear (:math:`^{13}\mathrm{C}`--:math:`^{1}\mathrm{H}`) dipolar
interactions, whereas carbonyl and quaternary carbons require the inclusion of
additional mechanisms such as chemical shift anisotropy (CSA), which arises
from the orientation dependence of the chemical shielding tensor and its
modulation by molecular rotational dynamics. Likewise, relaxation of
:math:`^{15}\mathrm{N}` nuclei generally requires a combined treatment of
dipolar and CSA mechanisms :cite:`allardNMRRelaxationMechanisms1997`.

Although these mechanisms are not currently implemented, extending
NMRDfromMD to support multi-nuclear relaxation analysis represents a natural
direction for future development.

Out-of-scope relaxation mechanisms
----------------------------------

While NMRDfromMD is designed to describe relaxation processes that can be
directly related to molecular motions extracted from classical MD trajectories,
some relaxation mechanisms require additional physical information beyond the
current framework. These mechanisms involve interactions that are not fully
captured by standard nuclear dipole-dipole relaxation models, such as
fluctuations of electric field gradients (EFG) or electron spin dynamics. The
following sections describe relaxation processes that are currently outside the
scope of NMRDfromMD and would require significant theoretical and
methodological extensions.

Quadrupolar relaxation
~~~~~~~~~~~~~~~~~~~~~~

NMRDfromMD is currently not intended for systems where quadrupolar
interactions represent the dominant relaxation mechanism. Nuclei with spin
quantum numbers :math:`I \geq 1`, such as :math:`^{2} \text{H}`,
:math:`^{14} \text{N}`, or other quadrupolar nuclei, require a different
theoretical treatment because quadrupolar relaxation originates from
fluctuations of the EFG tensor at the nucleus,
rather than from magnetic relaxation mechanisms such as dipole-dipole
interactions or CSA
:cite:`carofAccurateQuadrupolarNMR2014a, chubakQuadrupolar23NaNMR2023, vidalInitioMolecularDynamics2026`.

Predicting quadrupolar relaxation from MD trajectories requires the
calculation of the time-dependent EFG tensor along the trajectory, followed
by the evaluation of its autocorrelation function and corresponding spectral
densities. The EFG can be obtained using quantum-mechanical calculations or suitable
molecular models, for example classical force fields combined with an
appropriate description of electronic response, such as Sternheimer
antishielding corrections :cite:`sternheimerShieldingAntishieldingEffects1966`. These developments
introduce additional methodological and computational challenges and are
therefore outside the current scope of NMRDfromMD.

Paramagnetic relaxation
~~~~~~~~~~~~~~~~~~~~~~~

NMRDfromMD is currently not designed for systems in which paramagnetic
centers make a significant contribution to nuclear relaxation. In the
presence of unpaired electron spins, such as in systems containing
transition-metal ions or stable radicals, additional relaxation mechanisms
arise from electron--nuclear magnetic interactions. Unlike dipole--dipole or
CSA relaxation in diamagnetic systems, these mechanisms
depend not only on molecular dynamics but also on the dynamics of the
electron spin and its coupling to the surrounding nuclei
:cite:`cloreTheoryPracticeApplications2009b`.

Although molecular dynamics simulations can accurately describe structural
fluctuations and electron--nucleus distance distributions, they do not
provide the electronic properties required for quantitative predictions of
paramagnetic relaxation, such as electron spin relaxation times, hyperfine
coupling parameters, or magnetic susceptibility tensors. These quantities are
typically obtained from quantum-mechanical calculations, EPR measurements, or
dedicated spin dynamics approaches :cite:`westlundParamagneticEnhancedProton1993`.
Supporting paramagnetic relaxation would therefore require a specific theoretical
framework and is outside the current scope of NMRDfromMD.
