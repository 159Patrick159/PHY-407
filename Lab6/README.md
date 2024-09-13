# Lab 6: Molecular Dynamics
This sub-directory contains the python scripts and written report for Lab 6 for Computational Physics. The breakdown of the lab is the following:

## Q1 - Molecular Trajectories Simulation
Implement Verlet method to simulate system of two particles in cartestian coordinates subject to the Lennard-Jones potential. We perform three molecular simulations for different starting points for each particle. Lastly we verify energy conservation in our simulation by observing the KE and PE of each particle.

<p align='center'>
    <img src="Figures/Q1.png" title="Molecular simulation of particle trjactories subject to Lennard-Jones potential" height="80%" width="80%">
</p>

## Q2 - Molecular Dynamics Simulation
Building on our previous results we simulate the dynamic interactions of 16 particles that are being subejct to the Lennard-Jones potential. We observe there to be major energy fluctuaitons in the simulation which could be due to the implementation of the Verlet method.

<p align='center'>
    <img src="Figures/Q2a.png" title="Molecular simulation of dynamic interactions between 16 particles." height="80%" width="80%">
</p>