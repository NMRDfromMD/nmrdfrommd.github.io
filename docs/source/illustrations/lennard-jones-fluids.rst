.. include:: ../additional/links.rst
.. _lennard-jones-label:

Lennard-Jones fluid
===================

Here, NMR relaxation rates are measured from a Lennard-Jones fluid, and compared
to the results from Ref. :cite:`grivetNMRRelaxationParameters2005`.

MD system
---------

.. image:: lennard-jones-fluids/lj-dark.png
    :class: only-dark
    :alt: LJ fluid simulated with LAMMPS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

.. image:: lennard-jones-fluids/lj-light.png
    :class: only-light
    :alt: LJ fluid simulated with LAMMPS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

The system consists of 16,000 particles interacting via the classical 
Lennard-Jones (LJ) 12-6 potential and simulated using
LAMMPS :cite:`thompsonLAMMPSFlexibleSimulation2022`.
Each particle has a mass 
:math:`m = 1\,\text{g/mol}`, and LJ parameters 
:math:`\sigma = 3\,\text{Å}` and :math:`\epsilon = 0.1\,\text{kcal/mol}`.
All reduced parameters were taken to match the study by Grivet 
:cite:`grivetNMRRelaxationParameters2005`.
A cutoff of :math:`4 \sigma` was used for the LJ interactions. The 
simulation box has a volume of :math:`(26.9~\sigma)^3`, chosen to match 
the reduced density of :math:`\rho^* = 0.84`.
Production runs were performed in the microcanonical (NVE) ensemble, 
during which 10,000 timesteps were executed, equivalent to 50 times 
the reference time :math:`\sqrt{m \sigma^2/\epsilon}`. Configurations 
were recorded every 10 timesteps. A timestep of 
:math:`0.005\,\sqrt{m \sigma^2/\epsilon}` was used.
The imposed temperatures ranged from :math:`T = 30` to 
:math:`160\,\text{K}`, corresponding to reduced temperatures from 
:math:`T^* = 0.8` to :math:`3.0`.

All LAMMPS input scripts and analysis scripts written in Python are provided
on GitHub; see |dataset-LJ-fluid|.

Results
-------

The correlation function :math:`G^{(0)}` was first extracted for two temperatures, :math:`T = 50`
and :math:`140\,\text{K}`, and compared with the correlation functions reported by Grivet :cite:`grivetNMRRelaxationParameters2005`.
Our results show an excellent agreement with the results from Grivet, thus validating the
NMR formalism used here as well as the LJ system and parameters. 

.. image:: ../figures/illustrations/lennard-jones-fluid/G_correlation-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/illustrations/lennard-jones-fluid/G_correlation-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: Correlation function :math:`G^{(0)}` as extracted from the LJ fluid simulation
    for two different temperatures, and compared with the data from Grivet :cite:`grivetNMRRelaxationParameters2005` (open symbols).

The NMR relaxation rates :math:`R_1`
and :math:`R_2` was also extracted for all the temperatures, at
a frequency :math:`f_0 = 150\,\text{GHz}`. Our results
show a good agreement with the data from Grivet :cite:`grivetNMRRelaxationParameters2005`.

.. image:: ../figures/illustrations/lennard-jones-fluid/R1_spectra-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/illustrations/lennard-jones-fluid/R1_spectra-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: NMR relaxation rates :math:`R_1`
    and :math:`R_2` at
    a frequency :math:`f_0 = 150\,\text{GHz}`. 
    The data from Grivet :cite:`grivetNMRRelaxationParameters2005` are shown with open symbols.