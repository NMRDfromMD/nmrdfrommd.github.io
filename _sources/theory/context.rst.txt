Context
=======

A key ingredient for an accurate description of :math:`^1\mathrm{H}`
nuclear spin relaxation in soft matter systems is a realistic
representation of the stochastic rotational and translational motions
of molecules. Classical molecular dynamics (MD) simulations naturally
provide access to these dynamics and have therefore become the
reference approach for studying NMR relaxation.
Initially developed for simple fluids such as Lennard-Jones liquids,
MD has been used to characterize dipolar relaxation mechanisms
:cite:`odeliusIntermolecularDipoleDipoleRelaxation1993, grivetNMRRelaxationParameters2005`.
It has since been extended to increasingly realistic systems,
including water and other molecular liquids
:cite:`lippensT1RelaxationTime1993, calero1HNuclearSpin2015, singerMolecularDynamicsSimulations2017, singerNMRSpinrotationRelaxation2018, philipsProtonNMRRelaxation2019, singerElucidating1NMR2020`,
confined fluids in nanoporous materials
:cite:`khudozhitkovDynamicsPropenePropane2020, gravelleNMRInvestigationWater2023`,
and more complex soft matter systems such as polymers, lipid
membranes, proteins, and glass-forming liquids
:cite:`pastorLipidBilayersNMR2002, klaudaCollectiveNoncollectiveModels2008, singerElucidating1NMR2020, becherMolecularDynamicsSimulations2021, stenstromHowDoesIt2022`.

Beyond classical MD, other simulation techniques have also been
employed. Ab initio molecular dynamics has been used when electronic
structure effects are important, for example to study quadrupolar
relaxation
:cite:`calero1HNuclearSpin2015, philipsQuadrupolarNMRRelaxation2020, chubakNMRRelaxationRates2021`.
Monte Carlo simulations have also been explored
:cite:`friesMonteCarloCalculation1983`, although care must be taken
when extracting time-dependent correlation functions from
non-dynamical trajectories :cite:`huitemaCanMonteCarlo1999`.
More recently, coarse-grained models combined with structural
backmapping have been shown to reproduce NMR relaxation observables
:cite:`gravelleAssessingValidityNMR2023`.

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
