#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "src/ewald.c"
#include "src/constants.c"

#define TOL 1e-10

int main(void)
{
    const int n_particles = 500;
    const double box_size = 10.0;
    const int total_size = 3 * n_particles;

    double *pos = malloc(total_size * sizeof(double));
    double *pos_new = malloc(total_size * sizeof(double));
    double *charge = malloc(n_particles * sizeof(double));

    if (!pos || !pos_new || !charge)
        exit(EXIT_FAILURE);

    srand48(12345);

    /* ------------------ init system ------------------ */
    for (int i = 0; i < n_particles; i++)
    {
        pos[3 * i + 0] = drand48() * box_size;
        pos[3 * i + 1] = drand48() * box_size;
        pos[3 * i + 2] = drand48() * box_size;

        charge[i] = (i % 2 == 0) ? 1.0 : -1.0;
    }

    optimizeParameter(1e-4, box_size, charge, n_particles);

    /* ------------------ reference energy ------------------ */
    double E0 = ewd_reciprocal_space_coulomb_energy(
        pos, charge, n_particles, box_size);

    printf("\nE0 = %.15e\n", E0);

    /* ------------------ init S_k on OLD config ------------------ */
    init_Sk(pos, charge, n_particles, box_size);

    /* ------------------ multiple random moves ------------------ */
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

        /* --- brute-force delta --- */
        double E1 = ewd_reciprocal_space_coulomb_energy(
            pos_new, charge, n_particles, box_size);

        double dE_brute = E1 - E0;

        /* --- incremental delta --- */
        double dE_inc = delta_reciprocal_energy(
            i, pos_new, pos, charge, box_size);

        double err = fabs(dE_brute - dE_inc);

        printf("step %2d | dE_brute = %+ .12e | dE_inc = %+ .12e | err = %.3e\n",
               step, dE_brute, dE_inc, err);

        if (err > TOL)
        {
            printf("ERROR above tolerance!\n");
            return EXIT_FAILURE;
        }

        /* --- accept move, update everything --- */
        updateS_k(i, pos_new, pos, charge, box_size);
        memcpy(pos, pos_new, total_size * sizeof(double));
        E0 = E1;
    }

    printf("\nALL TESTS PASSED\n");

    free(pos);
    free(pos_new);
    free(charge);

    return 0;
}
