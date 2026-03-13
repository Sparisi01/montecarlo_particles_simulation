# Metropolis Monte Carlo Particle Simulations of Coupled Lennard-Jones and Coulomb Potential under Periodic Boundary Conditions

## Abstract

We present Monte Carlo simulations of a three-dimensional system of particles interacting via Lennard–Jones and Coulomb potentials under periodic boundary conditions. Equilibrium configurations are generated using the Metropolis algorithm within the canonical ensemble. Short-range interactions are treated using a cutoff scheme combined with a Verlet neighbor list, while long-range Coulomb interactions are computed via Ewald summation with optimized parameters. The implementation is validated against reference data for the pure Lennard–Jones fluid. We then investigate the effect of increasing ionic coupling strength on structural properties, analyzing energy and radial distribution functions. The results demonstrate the efficiency of the numerical approach and reveal a transition from liquid-like to strongly correlated ionic structures at high coupling strengths.

## Compilation and Usage

Compile and run the `main.c` file using the given `build_main.sh` file. Clang is required,

In the current form, all simulation setup has to be made inside the `main.c` file before compilation. In the following a list of variable that control the simulation:

- `double LAMBDA`: coupling constant between Coulomb and Lennard-Jones potential. If set to zero the Columb potential is turned off.
- `int lattice_type`: choose starting position lattice type (1 CC, 2 BCC, 4 FCC).
- `int n_cell_per_row`: number of lattice cell per row
- `double density`: particle numerical density
- `double temperature`: temperature used in the Metropolis Algorithm.
- `int N_thermalization_steps`: metropolis thermalization step.
- `int N_data_steps`: metropolis data step.

- `double VERLET_MAX_NEIGHTBOR_DISTANCE`: default value `3`.
- `double SKIN`: skinn radius in `VERLET_MAX_NEIGHTBOR_DISTANCE` units.
- `double ewald_error`: desired Ewald Summation statistical error.

At the end of the simulation `energy` and `radial_distribution` data are saved as `.csv` files in to the `output` folder. The `output` and `build` folders are automatically created by `build.sh` if not present.
