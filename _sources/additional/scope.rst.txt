Future development
==================

Currently, NMRDfromMD is primarily focused on the prediction and analysis
of :math:`^{1} \text{H}` nuclear magnetic relaxation dispersion (NMRD).
This focus is motivated by the relatively simple and well-established
relationship between proton relaxation rates and molecular motion, which is
often dominated by dipole-dipole interactions that can be directly described
through the autocorrelation functions of internuclear vectors extracted from
molecular dynamics trajectories. Extending the framework to other nuclei, such
as :math:`^{13} \text{C}`, :math:`^{15} \text{N}`, and :math:`^{19} \text{F}`,
is considerably more challenging because different relaxation mechanisms can
become dominant depending on the chemical environment. For example,
:math:`^{13} \text{C}` relaxation may be governed by heteronuclear
(:math:`^{13} \text{C}`-:math:`^{1} \text{H}`) dipolar interactions for
protonated carbons, whereas carbonyl or quaternary carbons require the
treatment of additional contributions such as chemical shift anisotropy (CSA)
and tensorial fluctuations. Similarly, relaxation analysis for nuclei such as
:math:`^{15} \text{N}` often requires the combined description of dipolar and
CSA mechanisms. These extensions therefore require additional theoretical
developments, including the calculation of anisotropic interaction tensors and
their time-dependent fluctuations from MD trajectories. Although not currently
implemented, expanding NMRDfromMD toward multi-nuclear relaxation analysis
represents a natural future direction for the package.

Scope limitations
-----------------

NMRDfromMD is currently not intended for systems where quadrupolar
interactions represent the dominant relaxation mechanism. Quadrupolar nuclei,
such as :math:`^{2} \text{H}`, :math:`^{14} \text{N}`, or other nuclei with
spin quantum numbers greater than :math:`1/2`, require a different theoretical
treatment because relaxation arises from fluctuations of the electric field
gradient (EFG) tensor rather than only from magnetic dipole-dipole
interactions or chemical shift anisotropy. Predicting quadrupolar relaxation
from MD trajectories would therefore require the calculation of the
time-dependent EFG tensor at the nucleus, either from suitable force-field
descriptions or from quantum-mechanical calculations, followed by the
evaluation of the corresponding spectral densities. These developments involve
additional methodological and computational challenges and are therefore
outside the current scope of NMRDfromMD.
