/**
 * @details
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>

#include "src/constants.c"
#include "src/checkpoints_handler.c"
#include "src/arrays_stat_operations.c"
#include "src/lennar_jones.c"
#include "src/periodic_boundaries.c"
#include "src/ewald.c"
#include "src/progress_bar.c"
#include "src/verlet_list.c"
#include "src/radial_distribution.c"

// Choose type of simulation
const double LAMBDA = 0; // For what lambda is see latex paper section "Units"
const int COULOMB_INTERACTION_ON = 1;

double pb_compute_total_energy(const double *pos_array,
                               const double *charge_array,
                               int n_particles,
                               int space_dim,
                               double box_size,
                               double epsilon,
                               double sigma)
{
    double total_energy = 0;
    total_energy += pb_total_lennar_jones_energy(pos_array, charge_array, n_particles, space_dim, box_size, epsilon, sigma);

    if (COULOMB_INTERACTION_ON)
    {
        if (space_dim != 3)
        {
            LOG_FATAL("Ewald Summation requires a space dimension of 3\n");
        }
        total_energy += LAMBDA * ewd_total_energy(pos_array, charge_array, n_particles, box_size);
    }

    return total_energy;
}

double pb_verlet_compute_total_energy(const double *pos_array,
                                      const double *charge_array,
                                      const IndexesList_t *verlet_list,
                                      int n_particles,
                                      int space_dim,
                                      double box_size,
                                      double epsilon,
                                      double sigma)
{
    double total_energy = 0;
    total_energy += pb_verlet_tot_lennar_jones_energy(pos_array, charge_array, verlet_list, n_particles, space_dim, box_size, epsilon, sigma);

    if (COULOMB_INTERACTION_ON)
    {
        if (space_dim != 3)
        {
            LOG_FATAL("Ewald Summation requires a space dimension of 3\n");
        }
        total_energy += LAMBDA * ewd_verlet_total_energy(pos_array, charge_array, verlet_list, n_particles, box_size);
    }

    return total_energy;
}

void init_system_random(double *pos_array,
                        double *charge_array,
                        int n_particles,
                        int space_dimension,
                        double box_size)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dimension; j++)
        {
            pos_array[c(i, j)] = drand48() * box_size; // Uniform position distribution inside the square box
        }

        switch (i % 2)
        {
        case 0:
            charge_array[i] = 1;
            break;

        default:
            charge_array[i] = -1;
            break;
        }
    }
}

int main(int argc, char const *argv[])
{
    FILE *complexity_file = fopen("./output/saves/Complexity/no_verlet_coul.csv", "w");

    int N_array[] = {100, 562, 1000, 1778, 3162, 5000, 7000, 10000, 15000, 25000};

    double density = 0.86;

    // In reduced unit keep those at 1
    double lennar_jones_epsilon = 1;
    double lennar_jones_sigma = 1;

    int space_dimension = 3; // 1D - 2D - 3D - ... - nD

    for (int q = 0; q < 10; q++)
    {

        int n_particles = N_array[q];
        double box_size = pow(n_particles / density, 1.0 / space_dimension);

        double space_step = 0.05;

        double VERLET_MAX_NEIGHTBOR_DISTANCE = 3 * lennar_jones_sigma;
        double SKIN = 0.5 * VERLET_MAX_NEIGHTBOR_DISTANCE;

        // Use the same r_c for verlet list and lennar jones
        LENNAR_JONES_CUT_OFF_IN_SIGMA_UNIT = VERLET_MAX_NEIGHTBOR_DISTANCE / lennar_jones_sigma;

        printf("N_particles : %d\n", n_particles);
        printf("Box size : %lf\n", box_size);

        int total_vel_pos_array_size = n_particles * space_dimension;

        double *pos_array = (double *)malloc(total_vel_pos_array_size * sizeof(double));

        double *charge_array = (double *)malloc(n_particles * sizeof(double));

        double energy = 0;

        init_system_random(pos_array, charge_array, n_particles, space_dimension, box_size);

        double *old_pos_array = (double *)malloc(total_vel_pos_array_size * sizeof(double));

        IndexesList_t *verlet_list = (IndexesList_t *)malloc(sizeof(IndexesList_t) * n_particles);

        // Build verlet list
        verlet_pb_build_list(pos_array, old_pos_array, verlet_list, n_particles, space_dimension, box_size, VERLET_MAX_NEIGHTBOR_DISTANCE, SKIN);

        printf("-------------------------------------\n");
        printf("VERLET_MAX_NEIGHTBOR_DISTANCE   : %.2f\n", VERLET_MAX_NEIGHTBOR_DISTANCE);
        printf("Max verlet count                : %d\n", get_max_verlet_count(verlet_list, n_particles));
        printf("-------------------------------------\n");

        // NOTE
        optimizeParameter(1, box_size, charge_array, n_particles);

        ewd_print_parameters();
        printf("-------------------------------------\n");

        REAL_CUTOFF = VERLET_MAX_NEIGHTBOR_DISTANCE;

        // Initialization of S(K) vector in the starting position
        ewd_init_S_K(pos_array, charge_array, n_particles, box_size);

        int N_step = 100;

        clock_t begin_time = clock();
        for (size_t i = 0; i < N_step; i++)
        {

            if (i % (N_step / 100) == 0)
            {

                // pb_verlet_compute_total_energy(pos_array, charge_array, verlet_list, n_particles, space_dimension, box_size, EPSILON, SIGMA);
                pb_compute_total_energy(pos_array, charge_array, n_particles, space_dimension, box_size, EPSILON, SIGMA);

                // Progress Bar
                print_progress(i, N_step, begin_time);
            }
        }

        clock_t end_time = clock();
        float seconds = (float)(end_time - begin_time) / CLOCKS_PER_SEC;
        // Clear terminal and print end simulation info
        printf("\r\033[2K");
        fprintf(complexity_file, "%d;%lf\n", n_particles, seconds / N_step);
        fflush(complexity_file);
        /* code */

        free(pos_array);
        free(charge_array);
        free(verlet_list);
        free(old_pos_array);
    }
}
