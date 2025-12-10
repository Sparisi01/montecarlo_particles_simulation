#ifndef CONSTANTS
#define CONSTANTS

#define SEED 42

#define N 500
#define BOX_SIZE 100
#define SPACE_DIM 2
#define N_METROPOLIS_STEPS 50000

#define PRINT_INTERVAL (N_METROPOLIS_STEPS / 100) // Update bar every 1%

#define PI 3.14159265359

#define c(i, j) (SPACE_DIM * i + j)

#endif