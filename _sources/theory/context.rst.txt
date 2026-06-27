
Context
=======

The measurement of NMR relaxation quantities allow for detailed studies of molecular motions
on time scales ranging from microseconds to minutes in systems as diverse as gases,
liquids, gels, polymers, adsorbed liquids, or solids
:cite:`goreNMRRelaxationWater1989, greiner-schmidSelfDiffusionCompressed1991`,
as well as proteins and other biological systems
:cite:`jacobsonProtonMagneticResonance1954, rorschachProteinDynamicsNMR1986`.

A key ingredient for an accurate description of nuclear spin relaxation of
:math:`^1 \text{H}` in soft matter systems is a realistic representation of
the stochastic rotational and translational motions of molecules. Classical
molecular dynamics (MD) simulations naturally provide access to these
dynamics, making them a widely used tool for studying NMR relaxation.

For simple fluids such as Lennard-Jones systems, MD has been used to
characterize dipolar relaxation mechanisms
:cite:`odeliusIntermolecularDipoleDipoleRelaxation1993, grivetNMRRelaxationParameters2005`. For more realistic molecular systems,
MD has been applied to water and small molecules
:cite:`lippensT1RelaxationTime1993, calero1HNuclearSpin2015, singerMolecularDynamicsSimulations2017, singerNMRSpinrotationRelaxation2018, philipsProtonNMRRelaxation2019, singerElucidatingNMRRelaxation2020`,
as well as to confined fluids in nanoporous materials
:cite:`khudozhitkovDynamicsPropenePropane2020, gravelleNMRInvestigationWater2023`. It has also been used for polymers,
lipid membranes, proteins, and glass-forming liquids such as glycerol
:cite:`becherMolecularDynamicsSimulations2021`.

Beyond classical MD, ab initio molecular dynamics has been employed to
compute NMR relaxation properties, particularly in cases where electronic
structure effects are important, such as quadrupolar relaxation
mechanisms :cite:`calero1HNuclearSpin2015,
philipsQuadrupolarNMRRelaxation2020,
chubakNMRRelaxationRates2021`. Monte Carlo simulations have also been
used :cite:`friesMonteCarloCalculation1983`, although care must be taken
when extracting time-dependent correlation functions from non-dynamical
trajectories :cite:`huitemaCanMonteCarlo1999`. Coarse-grained
models combined with structural backmapping have been shown to reproduce
NMR relaxation observables :cite:`gravelleAssessingValidityNMR2023`.

Despite the breadth of existing work, publicly available codes for computing NMR
relaxation from MD trajectories remain scarce. This limits reproducibility and
makes it difficult to apply established methods to new systems without significant
reimplementation effort. 

NMRDfromMD addresses this gap by providing an open-source, general-purpose code
for extracting NMR relaxation quantities directly from molecular dynamics
trajectories. It is designed to work with any MD engine capable of producing
standard trajectory formats, and covers isotropic liquids, polymer solutions,
and confined fluids. Numerical correctness and reproducibility are ensured
through a series of automated tests validated against well-established reference systems.

