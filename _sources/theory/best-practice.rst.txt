.. include:: ../additional/links.rst
.. _best-practice:

Best practices
==============

Accurate NMR relaxation calculations from molecular dynamics simulations require
careful attention to both the simulation protocol and the subsequent analysis.
Because relaxation rates depend on molecular structure and dynamics over a broad
range of timescales, they are sensitive to simulation parameters such as the force
field, trajectory length, sampling frequency, simulation box size, and analysis
settings. Here, some of the main factors that influence the accuracy of
NMR relaxation calculations from molecular dynamics are discussed, and practical
recommendations for obtaining reliable and reproducible results are provided.

Force field
-----------

The agreement between experiments and simulations is limited by the
quality of the chosen force field. While some force fields show excellent
agreement with NMR experimental data, for instance in simulations of water,
hydrocarbons, or polymer melts
:cite:`singerMolecularDynamicsSimulations2017,gravelleNMRInvestigationWater2023,gravelleAssessingValidityNMR2023`,
it is important to recognize that they are generally not parametrized
specifically to reproduce NMR relaxation observables. Instead, they are
typically optimized for selected thermodynamic, structural, and, in some
cases, dynamical properties
:cite:`mackerellEmpiricalForceFields2004,leachMolecularModellingPrinciples2001a`.
These target properties may include densities, heats of vaporization, phase
equilibria, solvation energies, radial distribution functions, or dynamical
observables. Since NMR relaxation rates are governed by time correlation
functions that reflect both equilibrium structure and molecular dynamics, the
suitability of a force field for relaxation studies depends on its ability to
capture the relevant motions on the corresponding timescales.

Impact of the water model
~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, the NMR relaxation properties of bulk water were calculated
from molecular dynamics simulations for three water models,
:math:`\text{TIP4P/2005}` :cite:`abascalGeneralPurposeModel2005`,
:math:`\text{SPC/E}` :cite:`berendsenMissingTermEffective1987`, and
:math:`\text{TIP3P}` :cite:`jorgensenComparisonSimplePotential1983`.

Correlation functions extracted for all three models at
:math:`T = 300\,\text{K}` highlight the differences between them, with
:math:`\text{TIP3P}` showing much faster decorrelation compared to
:math:`\text{SPC/E}`, which itself decorrelates slightly faster than
:math:`\text{TIP4P-2005}` (:ref:`Fig. 1 <fig:water-ff>`, panel A).
This ordering of the molecular dynamics timescales is consistent with the
relative viscosities reported for these models :cite:`gonzalezShearViscosityRigid2010` at :math:`T = 298\,\text{K}`:
:math:`0.321 \, \text{mPa s}` for :math:`\text{TIP3P}`,
:math:`0.729 \, \text{mPa s}` for :math:`\text{SPC/E}`, and
:math:`0.855 \, \text{mPa s}` for :math:`\text{TIP4P/2005}`.
Among the models considered, the :math:`\text{TIP4P/2005}` model has a viscosity closest to
the experimental value of :math:`0.896 \, \text{mPa s}`
:cite:`harrisTemperatureVolumeDependence2004`.

The NMR relaxation rate, :math:`R_1`, was extracted as a function of
temperature. The :math:`\text{TIP4P/2005}` model shows excellent agreement
with experimental measurements reported by Krynicki
:cite:`krynickiProtonSpinlatticeRelaxation1966` and Hindman et al.
:cite:`hindmanRelaxationProcessesWater2003`. By contrast, both
:math:`\text{SPC/E}` and :math:`\text{TIP3P}` underestimate the relaxation rate,
consistent with their larger deviations in viscosity and with previous
observations by Calero et al. :cite:`calero1HNuclearSpin2015`.

.. _fig:water-ff:

.. container:: figure

    .. image:: best-practice/water-ff-dm.png
        :class: only-dark
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. image:: best-practice/water-ff.png
        :class: only-light
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. container:: figurelegend

        Figure 1: A) Correlation functions extracted from molecular dynamics
        simulations for three water models:
        :math:`\text{TIP4P/2005}` (blue disks), :math:`\text{SPC/E}` (cyan squares),
        and :math:`\text{TIP3P}` (green pentagons) at
        :math:`T = 300\,\text{K}`.
        B) NMR relaxation rate, :math:`R_1`, as a function of the scaled inverse
        temperature, :math:`1000 [\text{K}]/T`, for
        bulk water obtained from molecular dynamics simulations using the
        :math:`\text{TIP4P/2005}`, :math:`\text{SPC/E}`, and
        :math:`\text{TIP3P}` models. Simulation results are compared with
        experimental measurements reported by Krynicki
        :cite:`krynickiProtonSpinlatticeRelaxation1966` and Hindman et al.
        :cite:`hindmanRelaxationProcessesWater2003`.

The water example illustrates a more general principle: NMR relaxation
provides a demanding test of molecular dynamics force fields because it is
sensitive to both equilibrium structure and molecular motions over a wide
range of timescales. Agreement with thermodynamic or structural properties
alone does not guarantee accurate relaxation rates. Consequently,
NMR relaxation calculations from MD simulations can be used not only to interpret
experimental measurements, but also to assess and compare force fields based on
their ability to reproduce dynamical observables.

Simulation protocol
-------------------

NMR relaxation calculations are sensitive to both thermodynamic and dynamical
properties. To ensure accurate simulations, the simulation protocol must be
carefully chosen alongside the force field discussed above. Important aspects
of the simulation protocol include the integration timestep, cutoff distances,
thermostat and barostat settings, and the equilibration procedure
:cite:`frenkelUnderstandingMolecularSimulation2002,allenComputerSimulationLiquids2017`.
An excessively large integration timestep introduces errors in the equations of
motion, insufficient sampling can bias calculated properties, and inappropriate
thermostat or barostat coupling parameters can artificially affect dynamical
properties. These effects can alter the structural
and dynamical properties of the system, which are directly reflected in the
time correlation functions used to calculate NMR relaxation rates and can
therefore lead to inaccurate relaxation predictions.

Cutoff
~~~~~~

The treatment of non-bonded interactions is an important consideration for NMR
relaxation calculations, as these interactions influence both the equilibrium
structure and molecular dynamics. These properties influence the time
correlation functions underlying relaxation rates. To quantify this effect, the NMR
relaxation rate :math:`R_1` of bulk water was calculated for different
Lennard-Jones cutoff distances :math:`r_\text{LJ}` using the
:math:`\text{TIP4P/2005}` model.

The intermolecular characteristic time, :math:`\tau_\text{T}`, shows a clear dependence
on the cutoff distance. For the smallest cutoff considered
(:math:`r_\text{LJ} = 0.6\,\text{nm}`), :math:`\tau_\text{T}` increases by about
:math:`10\,\%` compared to the converged value obtained for the largest cutoff
distances (:ref:`Fig. 2 <fig:water-co>`, panel A). This reflects
the weakening of intermolecular interactions caused by the truncation of the Lennard-Jones
potential, which modifies the liquid structure and the molecular motions contributing
to the intermolecular correlation functions. For a cutoff of
:math:`1\,\text{nm}`, which is a commonly used value in molecular dynamics
simulations, the deviation is reduced to approximately :math:`2\,\%`.

The variation of :math:`\tau_\text{T}` affects the calculated
intermolecular contribution to the relaxation rate,
:math:`R_\text{1, T}`, which is overestimated for the smallest
cutoffs (:ref:`Fig. 2 <fig:water-co>`, panel B). These observations are consistent
with previous measurements, see for instance Ref. :cite:`gravelleNMRInvestigationWater2023`.

.. _fig:water-co:

.. container:: figure

    .. image:: best-practice/water-co-dm.png
        :class: only-dark
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. image:: best-practice/water-co.png
        :class: only-light
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. container:: figurelegend

        Figure 2: a) Inter-molecular characteristic time, :math:`\tau_\text{T}`,
        as a function of the LJ cutoff, :math:`r_\text{LJ}`. The dashed line is a
        guide to the eye, indicating the value :math:`\tau_\text{T} = 3.78 \, \text{ps}`
        obtained for the largest cutoffs.
        b) Inter-molecular NMR relaxation rate, :math:`R_\text{1, T}`,
        as a function of :math:`r_\text{LJ}`.  The dashed line
        is :math:`R_\text{1, T} = 0.125 \, \text{s}^{-1}`.

Integration timestep
~~~~~~~~~~~~~~~~~~~~

The integration timestep determines the numerical accuracy of the molecular
dynamics trajectory. If the timestep is too large, the equations of motion are
not accurately integrated, resulting in systematic errors in both structural
and dynamical properties
:cite:`izaguirreLongerTimeSteps1999`.
Since NMR relaxation rates are directly related to
molecular motions, these integration errors can propagate into the calculated
correlation functions and relaxation rates. The timestep should therefore be
chosen according to established best practices for the selected force field and
simulation conditions, and its adequacy should be verified through convergence
testing when high accuracy is required.

Thermostat
~~~~~~~~~~

The thermostat controls the temperature of the simulated system, but it can
also influence molecular dynamics if applied too aggressively
:cite:`noseMolecularDynamicsMethod1984, berendsenMolecularDynamicsCoupling1984, hooverCanonicalDynamicsEquilibrium1985, tuckermanLiouvilleoperatorDerivedMeasurepreserving2006`.
Strong coupling or inappropriate thermostat parameters may artificially damp or modify
translational and rotational motions, leading to biased time-correlation
functions and relaxation rates. For NMR relaxation calculations, it is
therefore important to employ a thermostat that preserves realistic dynamics
and to use coupling parameters that minimally perturb the natural motion of
the system.

Simulation length
-----------------

The minimum simulation duration required to accurately calculate the NMR relaxation rate (e.g., :math:`R_1`)
depends on the quantity of interest. To obtain a converged value of :math:`R_1` in the zero-frequency limit,
the trajectory must be long enough for the correlation function
:math:`G(t)` to fully decay to zero, meaning the simulation duration
must significantly exceed the longest correlation time :math:`\tau_c`
of the system. For frequency-dependent quantities :math:`R_1(f)`, the
lowest accessible frequency is bounded by :math:`f_\text{min} \sim 1/T_\text{sim}`,
where :math:`T_\text{sim}` is the total simulation
duration. Frequencies below :math:`f_\text{min}` cannot be probed
regardless of the trajectory sampling interval. In practice, convergence
can be verified by comparing results obtained from simulations of
increasing duration.

Box size
--------

NMR relaxation measurements are sensitive to the 
finite-size effects that can occur with small simulation boxes :cite:`grivetNMRRelaxationParameters2005`.

As an illustration, the NMR relaxation rate :math:`R_1`
was measured for water with different number of molecules
:math:`N \in [100,\,10000]`, which correspond to equilibrium
box of lateral sizes :math:`L \in [1.4,\,6.7]\,\text{nm}`.
Our results show that the inter-molecular
relaxation rate :math:`R_1^\text{inter}` is sensitive to the 
box size even for the largest boxes considered here.
With small box size, the tail of :math:`G_\text{inter}`, 
which decreases as :math:`G_\text{inter} \sim t^{-3/2}`, is cutoff
which lead to an error on :math:`R_1^\text{inter}`.
Note that :math:`R_1^\text{intra}`, which is the dominant contribution to 
:math:`R_1` for water at ambient temperature, is barely affected
by the box size and therefore the resulting error induced on the 
total relaxation rate :math:`R_1` remains small for :math:`N > 1000`.

.. _fig:water-N:

.. container:: figure

    .. image:: best-practice/water-N-dm.png
        :class: only-dark
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. image:: best-practice/water-N.png
        :class: only-light
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. container:: figurelegend

        Figure 3: A) Inter-molecular NMR relaxation rate :math:`R_\text{1, T}` as a function of the number of molecules :math:`N`
        for a bulk water system. For the smallest systems, results were averaged
        from up to 10 independent simulations and the error bar is calculated from
        the standard deviation.
        The dashed line is a guide to the eye, indicating the value obtained
        for the largest value of :math:`N`.
        b) Inter-molecular correlation function :math:`G_{ij, \text{T}}`
        for two different numbers of molecules, :math:`N = 158` (cyan squares) and :math:`N = 10000` (green disks).

Trajectory output frequency
---------------------------

The trajectory output frequency sets the temporal resolution of the
analysis and determines the shortest correlation times that can be
resolved. The sampling interval :math:`\Delta t` must be significantly
smaller than the shortest relevant correlation time in the system,
otherwise fast molecular motions are not captured and the correlation
function :math:`G(t)` is undersampled. This leads to an overestimation
of the characteristic times :math:`\tau`, and consequently to errors
in the computed relaxation rates. If the characteristic correlation
times of the system are not known a priori, the appropriate
:math:`\Delta t` should be identified from convergence testing.
Note that a small sampling interval increases the size of the
trajectory files and the computational cost of the analysis.

As an illustration, the NMR relaxation time :math:`T_1` of bulk water
was measured for sampling intervals ranging from
:math:`\Delta t = 0.02\,\text{ps}` to :math:`5\,\text{ps}`. Using a
sampling interval larger than approximately :math:`\Delta t =
0.5\,\text{ps}` leads to a significant underestimation of :math:`T_1`,
accompanied by an overestimation of the inter-molecular characteristic
time :math:`\tau_\text{inter}`. Both effects are consistent with
insufficient sampling of fast rotational and translational motions of
the water molecules.

.. _fig:water-dump:

.. container:: figure

    .. image:: best-practice/water-dump-dm.png
        :class: only-dark
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. image:: best-practice/water-dump.png
        :class: only-light
        :alt: NMR results obtained from the LAMMPS simulation of water

    .. container:: figurelegend

        Figure 4: A) Intermolecular NMR relaxation time :math:`T_\text{1, T}`
        as a function of the trajectory dumping frequency :math:`\Delta t`
        for a bulk water system at :math:`T = 300 \text{K}`.
        The dashed line show the value for :math:`T_1`
        for :math:`\Delta t \to 0`.
        B) Intramolecular NMR relaxation time :math:`T_\text{1, R}`
        as a function of :math:`\Delta t`.

Analysis parameters
-------------------

The accuracy of NMR relaxation calculations depends not only on the quality
of the molecular dynamics simulation, but also on the parameters used during
the analysis. While these parameters do not alter the underlying trajectory,
they can affect the statistical uncertainty of the results or determine which
physical contributions are included in the calculation. Their influence should
therefore be assessed through appropriate convergence tests.

Number of reference atoms
~~~~~~~~~~~~~~~~~~~~~~~~~

The parameter ``number_i`` controls how many reference atoms are randomly
sampled during the calculation. Because the selection is stochastic,
results will vary slightly between runs when ``number_i > 0``. The
statistical uncertainty decreases as ``number_i`` increases, and
setting ``number_i = 0`` includes all eligible atoms, providing the
most accurate result at the highest computational cost. In practice,
convergence should be verified by repeating the calculation with
increasing values of ``number_i`` until the relaxation rates stabilize.

Cross-species interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``neighbor_group`` is not specified, intermolecular contributions
are computed only between atoms belonging to the same chemical species.
In a mixture such as polymer--water, this means that water--polymer
cross-interactions are excluded from the relaxation calculation. If
cross-species contributions are expected to be significant, the
appropriate ``neighbor_group`` must be set explicitly. Neglecting
these contributions may lead to an underestimation of the total
intermolecular relaxation rate.


.. list-table:: Summary of the main factors affecting the accuracy of NMR relaxation calculations from molecular dynamics.
   :header-rows: 1
   :widths: 25 35 40

   * - Parameter
     - If not properly chosen
     - Consequence

   * - Force field
     - Incorrect description of molecular structure or dynamics
     - Systematic deviations in relaxation rates due to inaccurate structural properties, molecular motions, or correlation times.

   * - Simulation protocol
     - Inappropriate timestep, cutoff distances, thermostat/barostat settings, or equilibration procedure
     - Errors in structural and dynamical properties, leading to inaccurate correlation functions and relaxation rates.

   * - Simulation length
     - Trajectory duration insufficient for correlation functions to fully decay
     - Incomplete convergence of relaxation rates, especially in the low-frequency limit.

   * - Box size
     - Simulation box too small, causing finite-size effects
     - Truncation of long-time intermolecular correlations and underestimation of intermolecular relaxation contributions.

   * - Trajectory output frequency
     - Sampling interval too large compared with relevant correlation times
     - Fast molecular motions are undersampled, causing errors in correlation functions and relaxation rates.

   * - Number of reference atoms
     - Too few atoms sampled during the calculation
     - Increased statistical uncertainty and reduced precision of the calculated relaxation rates.

   * - Cross-species interactions
     - Interactions between different chemical species are not included
     - Underestimation of intermolecular relaxation contributions in mixtures.
    