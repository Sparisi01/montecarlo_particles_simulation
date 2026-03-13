/**
 * @file test_ewald.c
 * @brief Test suite for Ewald summation implementation.
 *
 * This file contains a collection of numerical tests for validating
 * the correctness and stability of the Ewald summation algorithm
 * implemented in "src/ewald.c".
 *
 * Implemented tests:
 * - Convergence of total energy with respect to alpha parameter
 * - Energy invariance under spatial translation
 * - Consistency between brute-force and incremental long-range energy updates
 *
 * The tests validate:
 * - Numerical stability
 * - Translational invariance
 * - Correctness of delta-energy updates
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../src/ewald.c"
#include "../src/logger.c"

#define TOL 1e-10

int delta_test()
{
    LOG_TEST("------DELTA TEST------");

    const int n_particles = 500;
    const double box_size = 1;
    const int total_size = 3 * n_particles;

    ALPHA = 6;
    REAL_CUTOFF = box_size / 2;
    REAL_RANGE = 0;
    RECIPROCAL_RANGE = 6;

    double *pos = malloc(total_size * sizeof(double));
    double *pos_new = malloc(total_size * sizeof(double));
    double *charge = malloc(n_particles * sizeof(double));

    for (int i = 0; i < n_particles; i++)
    {
        pos[3 * i + 0] = drand48() * box_size;
        pos[3 * i + 1] = drand48() * box_size;
        pos[3 * i + 2] = drand48() * box_size;

        charge[i] = (i % 2 == 0) ? 1.0 : -1.0;
    }

    double E0 = ewd_long_energy(
        pos, charge, n_particles, box_size);

    printf("\nE0 = %.15e\n", E0);

    ewd_init_S_K(pos, charge, n_particles, box_size);

    for (int step = 0; step < 20; step++)
    {
        memcpy(pos_new, pos, total_size * sizeof(double));

        int i = lrand48() % n_particles;

        double dx = (drand48() - 0.5) * 0.2;
        double dy = (drand48() - 0.5) * 0.2;
        double dz = (drand48() - 0.5) * 0.2;

        pos_new[3 * i + 0] += dx;
        pos_new[3 * i + 1] += dy;
        pos_new[3 * i + 2] += dz;

        double E1 = ewd_long_energy(
            pos_new, charge, n_particles, box_size);

        double dE_brute = E1 - E0;

        double dE_inc = ewd_delta_long_energy(
            i, pos_new, pos, charge, box_size);

        double err = fabs(dE_brute - dE_inc);

        printf("step %2d | dE_brute = %+.12e | dE_inc = %+.12e | err = %+.3e\n",
               step, dE_brute, dE_inc, err);

        if (err > TOL)
        {
            printf("ERROR above tolerance!\n");
            LOG_TEST_FAILED;
        }

        ewd_update_S_K(i, pos_new, pos, charge, box_size);
        memcpy(pos, pos_new, total_size * sizeof(double));
        E0 = E1;
    }

    free(pos);
    free(pos_new);
    free(charge);

    LOG_TEST_PASSED;

    return 0;
}

void convergence_to_same_total_energy()
{

    LOG_TEST("------CONVERGENCE TEST------");

    FILE *file = fopen("./output/convergence_test_file.csv", "w");

    const double TOLLERANCE = 1e-6;

    const int n_particles = 100;
    const double box_size = 1;
    const int total_size = 3 * n_particles;

    double *pos = malloc(total_size * sizeof(double));
    double *charge = malloc(n_particles * sizeof(double));

    srand48(12345);

    for (int i = 0; i < n_particles; i++)
    {
        pos[3 * i + 0] = drand48() * box_size;
        pos[3 * i + 1] = drand48() * box_size;
        pos[3 * i + 2] = drand48() * box_size;

        charge[i] = (i % 2 == 0) ? 1.0 : -1.0;
    }

    const double alpha_min = 2;
    const double alpha_max = 40;
    const double n_alpha = 30;
    const double d_alpha = (alpha_max - alpha_min) / n_alpha;
    double old_en = 0;

    for (size_t i = 2; i < n_alpha; i++)
    {
        ALPHA = alpha_min + i * d_alpha;

        REAL_CUTOFF = 4.0 / ALPHA;
        REAL_RANGE = ceil(REAL_CUTOFF / box_size);

        RECIPROCAL_RANGE = ceil(8.0 * ALPHA * box_size / (2 * PI));

        double self_en = ewd_self_energy(charge, n_particles) / n_particles;
        double long_en = ewd_long_energy(pos, charge, n_particles, box_size) / n_particles;
        double short_en = ewd_short_energy(pos, charge, n_particles, box_size) / n_particles;
        double tot_en = long_en + short_en - self_en;

        printf("Alpha: %.1f, Short: %.5E, Long: %.5E, Self: %.5E, Tot: %.10E | ", ALPHA, short_en, long_en, self_en, tot_en);
        printf("alpha*rc = %.3f, kc/(2a)=%.3f\n",
               ALPHA * REAL_CUTOFF,
               (2 * PI * RECIPROCAL_RANGE / box_size) / (2 * ALPHA));
        fprintf(file, "%.f;%.5E;%.5E;%.5E;%.10E\n", ALPHA, short_en, long_en, self_en, tot_en);

        if (i > 4 && fabs(old_en - tot_en) / tot_en > TOLLERANCE)
        {
            LOG_TEST_FAILED;
            return;
        }

        old_en = tot_en;
    }

    LOG_TEST_PASSED;

    return;
}

void translate(double *pos, int n_particles, double dx, double dy, double dz, double L)
{
    for (int i = 0; i < n_particles; i++)
    {
        pos[3 * i + 0] = fmod(pos[3 * i + 0] + dx + L, L);
        pos[3 * i + 1] = fmod(pos[3 * i + 1] + dy + L, L);
        pos[3 * i + 2] = fmod(pos[3 * i + 2] + dz + L, L);
    }
}

// Check energy invariance under space-translation
void test_translation()
{

    LOG_TEST("------TRANSLATION TEST------");

    EWD_OPTIMIZED = 1;
    const int n_particles = 500;
    const double box_size = 1;
    const int total_size = 3 * n_particles;

    ALPHA = 6;
    REAL_CUTOFF = box_size / 2;
    REAL_RANGE = 0;
    RECIPROCAL_RANGE = 8;

    double *pos = malloc(total_size * sizeof(double));
    double *charge = malloc(n_particles * sizeof(double));

    for (int i = 0; i < n_particles; i++)
    {
        pos[3 * i + 0] = drand48() * box_size;
        pos[3 * i + 1] = drand48() * box_size;
        pos[3 * i + 2] = drand48() * box_size;

        charge[i] = (i % 2 == 0) ? 1.0 : -1.0;
    }

    double E1 = ewd_total_energy(pos, charge, n_particles, box_size);

    translate(pos, n_particles, 0.123, 0.456, 0.321, box_size);

    double E2 = ewd_total_energy(pos, charge, n_particles, box_size);

    printf("ΔE = %.3e\n", fabs(E1 - E2));

    if (fabs(E1 - E2) < TOL)
    {
        LOG_TEST_PASSED;
    }
    else
    {
        LOG_TEST_FAILED;
    }
}

int main()
{
    convergence_to_same_total_energy();
    test_translation();
    delta_test();

    printf("\n----------------\n");
    printf("ALL TESTS PASSED");
    printf("\n----------------\n");

    return 0;
}