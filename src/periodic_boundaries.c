#ifndef PERIODIC_BOUNDARIES
#define STRUCTURPERIODIC_BOUNDARIESES

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "constants.c"

/**
 * @brief Apply periodic boundary conditions to a position.
 *
 * Wraps a coordinate `x` into the primary simulation box of size `L`
 * using periodic boundary conditions.
 *
 * The returned value is guaranteed to lie in the interval:
 *     0 ≤ x < L
 *
 * This function is intended for wrapping particle positions, not
 * displacement vectors. It correctly handles values of `x` that are
 * outside the box by an arbitrary number of box lengths, including
 * negative values.
 *
 * Mathematical form:
 *     x_wrapped = x - floor(x / L) * L
 *
 * @param x  Position coordinate to be wrapped
 * @param box_size  Size of the periodic box (must be > 0)
 * @return   Wrapped position in the range [0, L)
 */
// ANCHOR - succesfully tested with 1 particle and constant force along x axis
static inline double pb_wrap_position(double x, double box_size)
{
    return x - floor(x / box_size) * box_size;
}

/**
 * @brief Apply the minimum image convention on position component dx.
 *
 * Given the displacement position dx between two particles,
 * applies periodic boundary conditions so that dx lies
 * in the range [-L/2, L/2].
 *
 * @param dx Displacement along x
 * @param box_size  Simulation box size (assumed cubic, L > 0)
 * @return   Minimum-image
 */
// ANCHOR - succesfully tested with 2 particle displaced in the corner of the box
static inline double pb_minimum_image(double dx, double box_size)
{

    dx -= floor(dx / box_size + 0.5) * box_size;
    return dx;
}

static inline double pb_pair_potential(double qi, double qj, double r2)
{
    // Prevent numerical blow-up
    if (r2 < 1e-32)
    {
        r2 = 1e-32;
    }

    double inv_r2 = 1. / r2;
    double inv_r = sqrt(inv_r2);
    double inv_r4 = inv_r2 * inv_r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    double inv_r12 = inv_r6 * inv_r6;

    double sigma_6 = SIGMA * SIGMA * SIGMA * SIGMA * SIGMA * SIGMA;
    double sigma_12 = sigma_6 * sigma_6;

    double V_coulomb = K_COUL * qi * qj * inv_r;
    // double V_Lennar_Jones = 4.0 * EPSILON * (sigma_12 * inv_r12 - sigma_6 * inv_r6);
    double V_8 = inv_r4 * inv_r4;

    return V_coulomb + V_8;
}

/**
 * @brief Compute the interaction energy of a single particle.
 *
 * Computes the total pairwise interaction energy between particle `i`
 * and all other particles in the system, applying periodic boundary
 * conditions with the minimum image convention.
 *
 * For each particle `k ≠ i`, the displacement vector is wrapped so that
 * each component lies in the range [-L/2, L/2], and the squared distance
 * r² is used to evaluate the pair potential.
 *
 * @param i             Index of the reference particle
 * @param pos_array     Array of particle positions (flattened)
 * @param charge_array  Array of particle charges
 * @param n_particles   Total number of particles
 * @param space_dim     Spatial dimension of the system
 * @param box_size      Size of the (cubic) periodic simulation box
 *
 * @return Interaction energy associated with particle `i`
 *
 * @note If this function is summed over all particles, the total energy
 *       will be double-counted and must be divided by two.
 */
double pb_compute_one_particle_energy(int i, const double *pos_array, const double *charge_array, int n_particles, int space_dim, double box_size)
{
    double energy_i = 0.0;

    for (int k = 0; k < n_particles; k++)
    {
        // Avoid self interaction
        if (i == k)
            continue;

        double r2 = 0.0;

        for (int j = 0; j < space_dim; j++)
        {
            double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
            r2 += dx * dx;
        }

        // Apply cutoff
        if (r2 > box_size * box_size)
            continue;

        energy_i += pb_pair_potential(charge_array[i], charge_array[k], r2);
    }

    return energy_i;
}

double pb_compute_tot_energy()
{
    return 0;
}

void pb_metropolis_step(double *pos_array, const double *charge_array, double delta, double temperature, int n_particles, int space_dim, long *accepted_counter, double box_size)
{
    // Alocate an array of dj on the stack,
    // this is used to keep track of the different direction and make the code independent on space dimension
    double dj_array[space_dim];

    // Update one particle at the time, the order of updating is fixed.
    for (int i = 0; i < n_particles; i++)
    {

        double old_i_energy = pb_compute_one_particle_energy(i, pos_array, charge_array, n_particles, space_dim, box_size);

        for (int j = 0; j < space_dim; j++)
        {
            // Random step in j direction between -delta and + delta
            double dj = ((drand48() * 2) - 1) * delta;

            // NOTE: In periodic boundaries conditions there is no need to check for boundaries
            dj_array[j] = dj;
            pos_array[c(i, j)] += dj;
        }

        double new_i_energy = pb_compute_one_particle_energy(i, pos_array, charge_array, n_particles, space_dim, box_size);
        double dE = new_i_energy - old_i_energy;

        // METROPOLIS ACCEPTANCE AND UPDATE energy_array
        int accepted = 0;
        double alpha = fmin(1, exp(-dE / temperature));
        // printf("%lf\n", alpha);
        accepted = drand48() <= alpha;

        // If the step is not accepted cancel the position update for the particle i
        if (!accepted)
        {
            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] -= dj_array[j];
            }
        }
        else
        {
            (*accepted_counter)++;

            // Bring the particle back to the first cell
            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] = pb_wrap_position(pos_array[c(i, j)], box_size);
            }
        }
    }
}

double pb_metropolis_step_all_system(double old_energy, double *pos_array, const double *charge_array, double delta, double temperature, int n_particles, int space_dim, long *accepted_counter, double box_size)
{
    // Save old position configuration
    // NOTE - This memory gets never free
    static double *old_pos_array = NULL;
    size_t total_array_size = sizeof(double) * n_particles * space_dim;

    if (old_pos_array == NULL)
    {
        old_pos_array = (double *)malloc(total_array_size);
        if (old_pos_array == NULL)
        {
            fprintf(stderr, "Memory allocation failed for old_pos_array\n");
            return -1;
        }
    }

    // Copy positions from pos_array to old_pos_array
    memcpy(old_pos_array, pos_array, total_array_size);

    // Update system: try to move each particle
    for (int i = 0; i < n_particles; i++)
    {
        for (int j = 0; j < space_dim; j++)
        {
            // Random step in j direction between -delta and + delta
            double dj = ((drand48() * 2) - 1) * delta;

            // NOTE: In periodic boundaries conditions there is no need to check for boundaries
            pos_array[c(i, j)] += dj;
        }
    }

    // Compute the new energy
    double new_energy = pb_compute_total_energy();
    double dE = new_energy - old_energy;

    // Metropolis acceptance criterion
    int accepted = 0;
    double alpha = fmin(1, exp(-dE / temperature));
    accepted = drand48() <= alpha;

    // If the step is not accepted cancel the position update for the particle i
    if (!accepted)
    {
        // Retrieve old positions
        memcpy(pos_array, old_pos_array, total_array_size);
    }
    else
    {
        (*accepted_counter)++;

        // Bring the particle back to the first cell
        for (size_t i = 0; i < n_particles; i++)
        {
            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] = pb_wrap_position(pos_array[c(i, j)], box_size);
            }
        }
    }

    return new_energy;
}

#endif