.. include:: ../additional/links.rst
.. _lennard-jones-label:

Code validation on a simple fluid
=================================

Here, the formalism presented in
:ref:`mathematical-framework` and implemented in NMRDfromMD is used
to predict the :math:`^1\mathrm{H}`-NMR relaxation properties of a
simple Lennard--Jones (LJ) fluid. The results are compared with the
reference calculations of Grivet :cite:`grivetNMRRelaxationParameters2005`.
The system consists of 16,000 particles interacting through the
classical LJ (12-6) potential, at a reduced density
:math:`\rho^* = 0.84` and reduced temperatures :math:`T^* = 0.8`\ -\ :math:`3.0`
(corresponding to :math:`T = 30`\ -\ :math:`160\,\text{K}`), matching the conditions of
Ref. :cite:`grivetNMRRelaxationParameters2005`. Simulation details are given
in :ref:`simulation-methods`; the input files and analysis scripts are available
on GitHub, see |dataset-LJ-fluid|.

Benchmark results
-----------------

The NMR relaxation rates :math:`R_1` and :math:`R_2` were evaluated at
a fixed frequency of :math:`f_0 = 150\,\mathrm{GHz}` (0.07 in reduced
units) over the full temperature range. The results show that
:math:`R_2` decreases monotonically from
:math:`5.6\,\mathrm{ms}^{-1}` at :math:`T = 30\,\text{K}` to
:math:`2.3\,\mathrm{ms}^{-1}` at :math:`T = 160\,\text{K}`. In
contrast, :math:`R_1` exhibits a maximum: it increases from
:math:`1.6\,\mathrm{ms}^{-1}` at :math:`T = 30\,\text{K}` to
:math:`1.9\,\mathrm{ms}^{-1}` around :math:`T = 80\,\text{K}`, before
decreasing to :math:`1.7\,\mathrm{ms}^{-1}` at
:math:`T = 160\,\text{K}` (:ref:`Figure 1 <fig:nmr-relaxation>`, panel A).
The results are in good agreement with the reference calculations of
Grivet :cite:`grivetNMRRelaxationParameters2005`, with typical
differences of approximately :math:`5\text{--}7\,\%`.

.. _fig:nmr-relaxation:

.. container:: figure

    .. image:: lennard-jones-fluids/nmr-relaxation-rates-at-target-dm.png
        :class: only-dark
        :alt: NMR relaxation rate of a LJ fluid simulated with LAMMPS
    
    .. image:: lennard-jones-fluids/nmr-relaxation-rates-at-target.png
        :class: only-light
        :alt: NMR relaxation rate of a LJ fluid simulated with LAMMPS
    
    .. container:: figurelegend

        Figure 1: A) NMR relaxation rates :math:`R_1` (circles) and
        :math:`R_2` (squares) as a function of temperature :math:`T`,
        computed from a Lennard--Jones simulation at a reduced frequency
        of 0.07, corresponding to
        :math:`f_0 = \omega_0/(2\pi) = 150\,\text{GHz}`.
        Results are compared with the reference data from
        Ref. :cite:`grivetNMRRelaxationParameters2005` (gray symbols).
        B) Product :math:`\omega_0 \tau` as a function of temperature.
        The dashed line indicates the condition
        :math:`\omega_0 \tau = 0.62`, corresponding to the BPP
        relaxation maximum. The inset shows a snapshot of the
        simulation, with atoms represented as spheres.

The temperature dependence of the relaxation rates can be understood in
terms of the evolution of the molecular correlation time, :math:`\tau`.
For this LJ fluid, the single-correlation-time approximation provides a
good description of the spectral density
:cite:`grivetNMRRelaxationParameters2005`. In this framework, the
spectral density follows a Lorentzian form:

.. math::
    :label: eq_J_Lorentzian

    J(\omega) = \frac{2 \tau}{1 + \omega^2 \tau^2}.

As shown in panel B, the correlation time decreases with increasing
temperature due to faster molecular motion. As a consequence,
:math:`J(0)`, which from :eq:`eq_J_Lorentzian` satisfies
:math:`J(0)=2\tau`, also decreases. A reduction of :math:`J(0)` leads to a
reduction of :math:`R_2` (:eq:`eq_BPP_R2`) and
explains the observed variation :math:`R_2` with :math:`T`.

In contrast, :math:`R_1` exhibits a non-monotonic temperature
dependence, with a maximum near :math:`T = 80`--:math:`90\,\text{K}`.
This behaviour is explained by the BPP relaxation model
(:eq:`eq_BPP_R1`). Inserting :eq:`eq_J_Lorentzian` into the BPP
expression for :math:`R_1`, the maximum occurs when

.. math::

    \omega \tau \approx 0.62.

The product :math:`\omega_0 \tau`, with :math:`\omega_0 = 2 \pi f_0` and :math:`\tau`
obtained directly from the simulations, crosses this value
around :math:`T = 80`--:math:`90\,\text{K}`,
which coincides with the maximum observed in :math:`R_1`
(:ref:`Figure 1 <fig:nmr-relaxation>`, panel B).
