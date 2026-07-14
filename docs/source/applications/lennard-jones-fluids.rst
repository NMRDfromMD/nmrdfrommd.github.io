.. include:: ../additional/links.rst
.. _lennard-jones-label:

Code validation on a simple fluid
==================================
 
This example validates NMRDfromMD by reproducing the :math:`^1\mathrm{H}`-NMR
relaxation properties of a Lennard--Jones (LJ) fluid, using the reference
calculations of Grivet :cite:`grivetNMRRelaxationParameters2005` as ground
truth. Because the system consists of identical spherical particles
interacting through a single distance-dependent pair potential, it provides a
clean, well-characterized benchmark, free of the ambiguities of molecular
geometry or intra-molecular contributions.

System
------

The system consists of 16,000 particles interacting through the classical
Lennard--Jones (12-6) potential, at a reduced density :math:`\rho^* = 0.84`
and reduced temperatures :math:`T^* = 0.8`\ -\ :math:`3.0`
(:math:`T = 30`\ -\ :math:`160\,\text{K}`), matching the conditions of Grivet
:cite:`grivetNMRRelaxationParameters2005`. Full simulation details are given
in :ref:`simulation-methods`; input and analysis scripts are available on
GitHub, see |dataset-LJ-fluid|.

Benchmark results
-----------------

We evaluate the NMR relaxation rates :math:`R_1` and :math:`R_2` at a fixed
frequency of :math:`f_0 = 150\,\mathrm{GHz}` (0.07 in reduced units) across
the full temperature range, and compare them with the reference data of
Grivet :cite:`grivetNMRRelaxationParameters2005`. Results show relatively good agreement
between the two data sets, with typical differences of about :math:`5-7\,\%`,

.. image:: lennard-jones-fluids/nmr-relaxation-rates-at-target-dm.png
    :class: only-dark
    :alt: NMR relaxation rate of a LJ fluid simulated with LAMMPS
 
.. image:: lennard-jones-fluids/nmr-relaxation-rates-at-target.png
    :class: only-light
    :alt: NMR relaxation rate of a LJ fluid simulated with LAMMPS
 
.. container:: figurelegend

    Figure: A) NMR relaxation rates :math:`R_1` (disks) and :math:`R_2` (squares) computed from
    the Lennard--Jones simulations at a frequency of 0.07 (dimensionless), or
    :math:`f_0 = 150\,\text{GHz}`, compared with the reference data from
    Grivet :cite:`grivetNMRRelaxationParameters2005` (gray symbols). B)
    Product :math:`\omega_0 \tau` as a function of temperature. The dashed
    line marks :math:`\omega_0 \tau = 0.62`, the value expected at the BPP
    relaxation maximum.

The relaxation rates exhibit the temperature dependence expected from the
shortening of the molecular correlation time, :math:`\tau`, as the fluid
becomes more mobile. :math:`R_2` decreases monotonically from
:math:`5.6\,\mathrm{ms}^{-1}` at :math:`T = 30\,\text{K}` to
:math:`2.3\,\mathrm{ms}^{-1}` at :math:`T = 160\,\text{K}`, reflecting the reduction of
the zero-frequency spectral density :math:`J(0)`, which dominates :math:`R_2`
but does not contribute to :math:`R_1`. :math:`R_1` instead passes through a
maximum: it rises from :math:`1.6\,\mathrm{ms}^{-1}` at :math:`T = 30\,\text{K}` to a
peak of :math:`1.9\,\mathrm{ms}^{-1}` around :math:`T = 80\,\text{K}`, then decreases
to :math:`1.7\,\mathrm{ms}^{-1}` at :math:`T = 160\,\text{K}`. This non-monotonic
behaviour is the signature of the BPP relaxation maximum, expected when the
correlation time satisfies :math:`\omega_0 \tau \approx 0.62`. Computing
:math:`\tau` directly from the simulations confirms this: :math:`\omega_0 \tau`
decreases from :math:`1.74` at :math:`T = 30\,\text{K}` to :math:`0.37` at
:math:`T = 160\,\text{K}`, crossing :math:`0.62` near 
:math:`T = 80\,\text{K}` and :math:`T = 90\,\text{K}`,
coinciding with the observed :math:`R_1` peak. At low
temperature (where :math:`\tau` is the largest), spectral density is concentrated near
zero frequency and contributes weakly to :math:`R_1`. As :math:`\tau`
decreases, spectral density shifts toward the Larmor frequency
:math:`\omega_0`, increasing :math:`R_1`, until :math:`\tau` becomes short
enough that the fluid enters the fast-motion (extreme narrowing) regime and
:math:`R_1` decreases again.
