Scope
=====

Current capabilities
--------------------

NMRDfromMD is currently designed for the prediction and analysis of 1H nuclear
magnetic relaxation dispersion (NMRD) profiles from molecular dynamics
trajectories. The current implementation focuses on relaxation mechanisms
dominated by magnetic dipole-dipole interactions, which can be directly
described through the time-dependent fluctuations of internuclear vectors. In
particular, 1H relaxation can be calculated from the rotational autocorrelation
functions of relevant interatomic vectors extracted from MD trajectories,
followed by the calculation of the corresponding spectral densities.

The current framework is applicable to both isotropic and anisotropic molecular
systems. For isotropic systems, such as bulk liquids or freely tumbling
molecules, the relaxation analysis relies on the assumption of isotropic
rotational diffusion. For anisotropic systems, such as confined liquids or
molecules interacting with surfaces, the full orientation-dependent dynamics of
internuclear vectors are retained, allowing the effect of restricted or
heterogeneous molecular motions on the relaxation dispersion profile to be
investigated. This enables the study of systems where deviations from simple
isotropic rotational diffusion play a key role in determining the observed NMR
relaxation behaviour.

Future extensions
-----------------

Currently, NMRDfromMD is primarily focused on the prediction and analysis of 1H
nuclear magnetic relaxation dispersion (NMRD). This focus is motivated by the
relatively simple and well-established relationship between proton relaxation
rates and molecular motion, which is often dominated by dipole-dipole
interactions that can be directly described through the autocorrelation
functions of internuclear vectors extracted from molecular dynamics
trajectories. Extending the framework to other nuclei, such as 13C, 15N, and
19F, is considerably more challenging because different relaxation mechanisms
can become dominant depending on the chemical environment. For example, 13C
relaxation may be governed by heteronuclear (13C-1H) dipolar interactions for
protonated carbons, whereas carbonyl or quaternary carbons require the treatment
of additional contributions such as chemical shift anisotropy (CSA) and
tensorial fluctuations. Similarly, relaxation analysis for nuclei such as 15N
often requires the combined description of dipolar and CSA mechanisms. These
extensions therefore require additional theoretical developments, including the
calculation of anisotropic interaction tensors and their time-dependent
fluctuations from MD trajectories. Although not currently implemented, expanding
NMRDfromMD toward multi-nuclear relaxation analysis represents a natural future
direction for the package.

Out-of-scope relaxation mechanisms
----------------------------------

While NMRDfromMD is designed to describe relaxation processes that can be
directly related to molecular motions extracted from classical MD trajectories,
some relaxation mechanisms require additional physical information beyond the
current framework. These mechanisms involve interactions that are not fully
captured by standard nuclear dipole-dipole relaxation models, such as
fluctuations of electric field gradients or electron spin dynamics. The
following sections describe relaxation processes that are currently outside the
scope of NMRDfromMD and would require significant theoretical and methodological
extensions.

Quadrupolar relaxation
~~~~~~~~~~~~~~~~~~~~~~

NMRDfromMD is currently not intended for systems where quadrupolar interactions
represent the dominant relaxation mechanism. Quadrupolar nuclei, such as 2H,
14N, or other nuclei with spin quantum numbers greater than 1/2, require a
different theoretical treatment because relaxation arises from fluctuations of
the electric field gradient (EFG) tensor rather than only from magnetic dipole-
dipole interactions or chemical shift anisotropy. Predicting quadrupolar
relaxation from MD trajectories would therefore require the calculation of the
time-dependent EFG tensor at the nucleus, either from suitable force-field
descriptions or from quantum-mechanical calculations, followed by the evaluation
of the corresponding spectral densities. These developments involve additional
methodological and computational challenges and are therefore outside the
current scope of NMRDfromMD.

Paramagnetic relaxation
~~~~~~~~~~~~~~~~~~~~~~~

NMRDfromMD is also currently not designed for systems where paramagnetic centers
provide a significant contribution to nuclear relaxation. In the presence of
unpaired electron spins, such as in systems containing transition-metal ions or
other paramagnetic species, relaxation can be dominated by electron-nuclear
dipolar interactions and contact interactions. These mechanisms depend on the
dynamics of the electron spin, the electron-nuclear distance distribution, and
the electron spin relaxation properties, which are not directly accessible from
standard classical molecular dynamics trajectories. A quantitative description
of paramagnetic relaxation would therefore require additional information, such
as electron spin relaxation times, hyperfine coupling parameters, or magnetic
susceptibility tensors, often obtained from quantum-mechanical calculations or
specialized spin dynamics approaches. Such extensions represent a major increase
in theoretical complexity and are therefore outside the current scope of
NMRDfromMD.
