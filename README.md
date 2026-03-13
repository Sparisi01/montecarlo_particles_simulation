# Metropolis Monte Carlo Simulation of Lennard-Jones and Coulomb Systems

![Build](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)

Monte Carlo simulation of interacting particles under **periodic boundary conditions** using the **Metropolis algorithm** in the canonical ensemble (NVT).

For a full project report see [report.pdf](report.pdf).

The code supports systems interacting through:

- **Lennard–Jones potential** (short-range interactions)
- **Coulomb potential** computed with **Ewald summation**

To improve computational efficiency, short-range interactions are evaluated using **Verlet neighbor lists**, reducing the complexity from O(N²) to approximately O(N) per Monte Carlo step.

The implementation is written in **C** and designed for studying the equilibrium properties of dense particle systems such as **Lennard-Jones fluids, ionic liquids, and strongly coupled plasmas**.

---

## Table of Contents

- [Features](#features)
- [Physical Model](#physical-model)
- [Units](#units)
- [Compilation](#compilation)
- [Main Simulation Parameters](#main-simulation-parameters)
  - [Interaction parameters](#interaction-parameters)
  - [System configuration](#system-configuration)
  - [Thermodynamic parameters](#thermodynamic-parameters)
  - [Monte Carlo parameters](#monte-carlo-parameters)
  - [Simulation length](#simulation-length)
  - [Neighbor list parameters](#neighbor-list-parameters)
  - [Ewald summation](#ewald-summation)
- [Output](#output)
  - [Energy time series](#energy-time-series)
  - [Radial distribution functions](#radial-distribution-functions)
- [Validation](#validation)
  - [Performance Scaling](#performance-scaling)
  - [Lennard-Jones simulation against NIST reference data](#lennard-jones-simulation-against-nist-reference-data)
- [Possible Extensions](#possible-extensions)

---

## Features

- Metropolis Monte Carlo sampling (NVT ensemble)
- Lennard-Jones potential with cutoff
- Tail correction for truncated Lennard-Jones interactions
- Coulomb interaction via **Ewald summation**
- Periodic boundary conditions
- **Verlet neighbor lists** for efficient short-range interactions
- Radial distribution function $g(r)$
- Charge-resolved radial distributions $g_{++}(r)$, $g_{--}(r)$, $g_{+-}(r)$
- Lattice initialization (CC, BCC, FCC)
- Binary checkpoint system for restarting simulations
- CSV output for easy post-processing

[⬆ Back to top](#table-of-contents)

---

## Physical Model

Particles interact through the potential

$U(r) = U_{LJ}(r) + \lambda U_{Coulomb}(r)$

### Lennard-Jones potential

$U_{LJ}(r) = 4\epsilon [ (\sigma/r)^{12} − (\sigma/r)^6 ]$

### Coulomb potential

$U_{Coulomb}(r) = \frac{q_i q_j}{r}$

The Coulomb interaction is computed using **Ewald summation**, which decomposes the interaction into:

- real-space short-range contribution
- reciprocal-space contribution
- self-interaction correction

[⬆ Back to top](#table-of-contents)

---

## Units

Simulations are performed in **Lennard-Jones reduced units**:

| Quantity     | Base unit               | Conversion                             |
|--------------|-------------------------|----------------------------------------|
| Length       | $\sigma$                | $r = r^* \sigma$                       |
| Energy       | $\varepsilon$           | $\mathcal{U} = \mathcal{U}^* \varepsilon$ |
| Temperature  | $\varepsilon/k_B$       | $T = T^* \frac{\varepsilon}{k_B}$      |
| Charge       | $e$                     | $q = q^* e$                            |

By introducing the dimensionless coupling constant $\lambda$ defined as

$\lambda \equiv \frac{e^2}{4\pi\epsilon_0 \sigma \varepsilon}$

the total pair potential becomes

$u^*_{tot}(r^*) = 4(\frac{1}{{r^*}^{12}} - \frac{1}{{r^*}^6}) + \lambda \frac{q^*_iq^*_j}{r^*}$

[⬆ Back to top](#table-of-contents)

---

## Compilation

Compile the code using the provided build script:

```t
./build_main.sh
```

Requirements:

- `clang`
- standard C library
- Unix-like environment

The build script automatically creates the directories:

```t
build/
output/
```

[⬆ Back to top](#table-of-contents)

---

## Main Simulation Parameters

The most important parameters are located at the beginning of `main.c`.

### Interaction parameters

```c
double LAMBDA
```

Coupling constant controlling the strength of the Coulomb interaction.

```c
LAMBDA = 0
```

disables Coulomb interactions (pure Lennard-Jones system).

---

### System configuration

```c
int lattice_type
```

Initial particle lattice.

| Value | Lattice |
|------|------|
| 1 | Simple cubic |
| 2 | BCC |
| 4 | FCC |

```c
int n_cell_per_row
```

Number of lattice cells per dimension.

Total number of particles:

```c
int n_particles = pow(n_cell_per_row, 3) * lattice_type;
```

---

### Thermodynamic parameters

```c
double density
double temperature
```

- particle number density
- Monte Carlo temperature

[⬆ Back to top](#table-of-contents)

---

## Monte Carlo parameters

```c
double space_step
```

Maximum displacement attempted during a Metropolis move.

Typical acceptance rate:

```t
40% – 70%
```

[⬆ Back to top](#table-of-contents)

---

## Simulation length

```c
int N_thermalization_steps
int N_data_steps
```

- thermalization steps
- production steps used for statistics

![equilibrization](png/equilibrization.png)

[⬆ Back to top](#table-of-contents)

---

## Neighbor list parameters

```c
double VERLET_MAX_NEIGHTBOR_DISTANCE
double SKIN
```

These define the cutoff radius and skin distance used for the Verlet neighbor list.

The list is rebuilt when particles move more than the **skin distance**.

[⬆ Back to top](#table-of-contents)

---

## Ewald summation

```c
double ewald_error
```

Target accuracy used to optimize Ewald parameters.

[⬆ Back to top](#table-of-contents)

---

## Output

Simulation results are written to the `output/` directory.

### Energy time series

```t
output/energy.csv
```

Contains the total system energy at each Monte Carlo step.

---

### Radial distribution functions

```t
output/radial_distribution.csv
output/radial_distribution_differ.csv
output/radial_distribution_equal.csv
```

These files contain:

| File | Description |
|-----|-----|
| radial_distribution | total $g(r)$ |
| radial_distribution_differ | $g_{+-}(r)$ |
| radial_distribution_equal | $g_{++}(r)$ and $g_{--}(r)$ |

The format is:

```t
r ; g(r)
```

and can be easily processed using Python.

[⬆ Back to top](#table-of-contents)

---

## Validation

### Performance Scaling

The following benchmark shows the runtime per Monte Carlo step as a function of the number of particles.

![Scaling benchmark](png/complexity_styled.png)

The measurements were performed at fixed particle density while increasing the system size. Each data point represents the average runtime per Monte Carlo step over a production run.

The near-linear behavior confirms the expected performance improvement provided by the neighbor list algorithm.

$O(N^{3/2})$ is the expected optimized behavior for classical Ewald summation algorithm.

### Lennard-Jones simulation against NIST reference data

Results from simulations of a pure Lennard–Jones (LJ) system compared with the NIST reference data (<https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm>).

The table reports, in order:

- temperature \(T\)
- density \(\rho\)
- mean energy per particle \(\mathcal{U}/N\)
- reference value from NIST
- acceptance probability \(P_A\)
- step autocorrelation time \(\tau\)

Energy errors are estimated as the standard deviation of the energy divided by the square root of the effective number of uncorrelated samples:

$N_{\text{eff}} = \frac{N_{\text{steps}}}{2\tau}$

The uncertainties correspond to a **68% confidence level**.

| T | ρ | U/N | NIST U/N | P_A | τ |
|---|---|---|---|---|---|
| 0.85 | 0.001 | -0.009799 ± 1.1×10⁻⁵ | -0.01032 ± 2×10⁻⁵ | 99% | 0.5 |
| 0.85 | 0.007 | -0.07214 ± 3×10⁻⁵ | -0.07283 ± 1.3×10⁻⁴ | 99% | 0.5 |
| 0.85 | 0.776 | -5.511 ± 2×10⁻³ | -5.5121 ± 4×10⁻⁴ | 50% | 206 |
| 0.85 | 0.86 | -6.023 ± 3×10⁻³ | -6.0305 ± 2.3×10⁻⁴ | 42% | 289 |

In the following the euilibrium $g(r)$ for all 4 systems:

![Scaling benchmark](png/g_dual_plot.png)

[⬆ Back to top](#table-of-contents)

---

## Possible Extensions

Potential improvements include:

- pressure estimator
- automatic Metropolis step tuning
- parallelization (OpenMP / MPI)
- molecular dynamics integration

[⬆ Back to top](#table-of-contents)
