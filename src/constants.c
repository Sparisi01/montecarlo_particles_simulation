#ifndef CONSTANTS
#define CONSTANTS

#define SEED 42
#define PRINT_INTERVAL (N_METROPOLIS_STEPS / 100) // Update bar every 1%
#define c(i, j) (SPACE_DIM * i + j)

// Simulation Parameters
#define N 50
#define SPACE_DIM 2
#define BOX_SIZE 100
#define N_METROPOLIS_STEPS 400000
#define STEP_SIZE 0.05
#define TEMPERATURE 0.1

// Physical constants
#define EPSILON 1  // Depth of LJ well
#define SIGMA 0.05 // Zero-crossing distance
#define K_COUL 1.0 // Coulomb prefactor (1/(4*pi*eps0) if SI)
#define PI 3.14159265359

#endif