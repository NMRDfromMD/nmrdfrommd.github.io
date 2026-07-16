Scope
=====

Current capabilities
--------------------

NMRDfromMD is currently designed for the prediction and analysis of
:math:`^{1} \text{H}` nuclear magnetic relaxation dispersion (NMRD) profiles
from molecular dynamics trajectories. The current implementation focuses on
relaxation mechanisms that can be described from time-dependent fluctuations
of molecular properties extracted directly from MD trajectories, from both isotropic and anisotropic
molecular systems :cite:`amaro-estradaQuantifyingEffectSpatial2022, gravelleIntermittentMolecularMotion2025`.
For protons, these fluctuations are typically dominated by dipole--dipole interactions.
:math:`^{1} \text{H}` relaxation can be calculated from the autocorrelation functions
of interatomic vectors extracted from MD trajectories, followed by spectral densities calculations
:cite:`singerMolecularDynamicsSimulations2017, valiyaparambathuMolecularModesElucidate2024`.

Future extensions
-----------------

The dominant relaxation mechanisms depend strongly on the nuclear species and
chemical environment. Therefore, extending the framework to nuclei other than
:math:`^{1} \text{H}`, such as :math:`^{13}\mathrm{C}`,
:math:`^{15}\mathrm{N}`, and :math:`^{19}\mathrm{F}`, requires the
implementation of additional relaxation mechanisms that may become important
for specific chemical environments. For example, :math:`^{13}\mathrm{C}`
relaxation in protonated carbons is often governed by heteronuclear
(:math:`^{13}\mathrm{C}` - :math:`^{1}\mathrm{H}`) dipole--dipole
interactions, whereas carbonyl and quaternary carbons may require the
inclusion of chemical shift anisotropy (CSA), which arises from the
orientation dependence of the chemical shielding tensor and its modulation
by molecular rotational dynamics :cite:`allardNMRRelaxationMechanisms1997`. Likewise, relaxation of :math:`^{15}\mathrm{N}`
nuclei generally requires a combined treatment of dipolar and CSA mechanisms
:cite:`fushmanDirectMeasurementOf151998`.

Although these mechanisms are not currently implemented, extending
NMRDfromMD to support multi-nuclear relaxation analysis represents a natural
direction for future development.

Out-of-scope relaxation mechanisms
----------------------------------

While NMRDfromMD is designed to describe relaxation processes that can be
related to molecular motions extracted from classical MD trajectories, some
relaxation mechanisms require additional physical information not typically
available from standard classical MD simulations. The following sections
describe relaxation processes that are currently outside the scope of
NMRDfromMD and would require significant theoretical and methodological
extensions.

Quadrupolar relaxation
~~~~~~~~~~~~~~~~~~~~~~

Nuclei with spin quantum numbers :math:`I \geq 1`, such as
:math:`^{2} \text{H}` and :math:`^{14} \text{N}`, often experience relaxation
dominated by quadrupolar interactions. Quadrupolar relaxation originates from
fluctuations of the electric field gradient (EFG) tensor at the nucleus
:cite:`carofAccurateQuadrupolarNMR2014a, chubakQuadrupolar23NaNMR2023, vidalInitioMolecularDynamics2026`. 
Therefore, predicting quadrupolar relaxation from MD trajectories requires the
calculation of the time-dependent EFG tensor along the trajectory, followed by the 
evaluation of its autocorrelation function and corresponding spectral densities. The EFG can be
obtained using quantum-mechanical calculations or suitable molecular models,
for example classical force fields combined with an appropriate description
of electronic response, such as Sternheimer antishielding corrections
:cite:`sternheimerShieldingAntishieldingEffects1966`. These developments
introduce additional methodological requirements and are therefore outside the
current scope of NMRDfromMD.

Paramagnetic relaxation
~~~~~~~~~~~~~~~~~~~~~~~

In the presence of paramagnetic centers containing unpaired electrons, such as
transition-metal ions or stable radicals, additional relaxation mechanisms
arise from electron--nuclear magnetic interactions. These mechanisms depend
not only on molecular dynamics but also on electron spin dynamics and its
coupling to surrounding nuclei
:cite:`cloreTheoryPracticeApplications2009b`. Classical MD simulations can
describe structural fluctuations and electron--nucleus distance distributions,
but they do not typically provide the electronic properties required for
quantitative predictions of paramagnetic relaxation, such as electron spin
relaxation times, hyperfine coupling parameters, or magnetic susceptibility
tensors. These quantities can be obtained from quantum-mechanical
calculations or dedicated spin dynamics approaches
:cite:`westlundParamagneticEnhancedProton1993`. Supporting paramagnetic
relaxation would therefore require additional theoretical developments beyond
the current MD-based relaxation framework.
