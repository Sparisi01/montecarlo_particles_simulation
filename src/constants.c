#ifndef CONSTANTS
#define CONSTANTS

#define SEED 42

#define c(i, j) (SPACE_DIM * i + j)

// Simulation Parameters
#define N 100
#define SPACE_DIM 3
#define BOX_SIZE 10
#define N_METROPOLIS_STEPS 100000
#define STEP_SIZE 0.1
#define TEMPERATURE 100

// Physical constants
#define EPSILON 1. // Depth of LJ well
#define SIGMA 1.   // Zero-crossing distance
#define K_COUL 1.  // Coulomb prefactor (1/(4*pi*eps0) if SI)
#define PI 3.14159265359
#define SQRT_PI 1.772453851

#endif