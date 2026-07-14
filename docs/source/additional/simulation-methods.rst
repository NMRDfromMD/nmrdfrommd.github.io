.. _simulation-methods:
 
Simulation methods
===================
 
This page collects the molecular dynamics (MD) simulation details for the
systems presented throughout this documentation. Each application or
validation page links back here rather than repeating the full setup, so
that simulation parameters are documented in a single, consistent location.

Lennard-Jones fluid
--------------------

.. image:: ../applications/lennard-jones-fluids/lj-dark.png
    :class: only-dark
    :alt: LJ fluid simulated with LAMMPS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right
 
.. image:: ../applications/lennard-jones-fluids/lj-light.png
    :class: only-light
    :alt: LJ fluid simulated with LAMMPS - Dipolar NMR relaxation time calculation
    :width: 250
    :align: right

The simulated system contains 16,000 particles interacting through the
classical Lennard--Jones (12-6) potential and was simulated using
LAMMPS :cite:`thompsonLAMMPSFlexibleSimulation2022`. Each particle has a
mass :math:`m = 1\,\mathrm{g/mol}` together with LJ parameters
:math:`\sigma = 3\,\text{Å}` and :math:`\epsilon = 0.1\,\mathrm{kcal/mol}`.
All reduced simulation parameters were chosen to reproduce the study of
Grivet :cite:`grivetNMRRelaxationParameters2005`. In particular, the
interaction cutoff was set to :math:`4 \sigma`, while the cubic
simulation box had a side length of :math:`26.9 \sigma`, corresponding to
the reduced density :math:`\rho^*=0.84`. Production runs were performed in the
microcanonical (NVE) ensemble, during which 10,000 timesteps were executed,
equivalent to 50 times the reference time :math:`\sqrt{m \sigma^2/\epsilon}`.
Configurations were recorded every 10 timesteps. A timestep of 
:math:`0.005\,\sqrt{m \sigma^2/\epsilon}` was used.
The imposed temperatures ranged from :math:`T = 30` to 
:math:`160\,\text{K}`, corresponding to reduced temperatures from 
:math:`T^* = 0.8` to :math:`3.0`.

.. admonition:: Reduced Lennard--Jones units
    :class: non-title-info

    Lennard-Jones simulations are commonly expressed in reduced units,
    where the particle mass :math:`m`, the characteristic length
    :math:`\sigma`, and the interaction energy :math:`\epsilon` define
    the natural scales of the system. Using reduced variables allows
    simulations with different physical parameters to be compared
    directly.