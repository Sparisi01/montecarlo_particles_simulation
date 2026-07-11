# Metropolis Monte Carlo Particle Simulations of Coupled Lennard-Jones and Coulomb Potential under Periodic Boundary Conditions

![Build](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)

Monte Carlo simulation of interacting particles under **periodic boundary conditions** using the **Metropolis algorithm** in the canonical ensemble (NVT).

For a full project report see [report.pdf](report.pdf).

The code supports systems interacting through:

- **Lennard–Jones potential** (short-range interactions)
- **Coulomb potential** computed with **Ewald summation**

To improve computational efficiency, short-range interactions are evaluated using **Verlet neighbor lists**, reducing the complexity from O(N²) to approximately O(N) per Monte Carlo step.

The implementation is written in **C** and designed for studying the equilibrium properties of dense particle systems such as **Lennard-Jones fluids, ionic liquids, and strongly coupled plasmas**.

The implementation is succesfully tested again the standard NIST (SRD 173) Lennard-Jones fluid simulations database.
A study of the equilibrium properties varying the LJ and Coulomb coupling constant is performed.
