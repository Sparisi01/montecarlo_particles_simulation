#ifndef LENNAR_JONES_C
#define LENNAR_JONES_C

#include <stdlib.h>
#include "verlet_list.c"
#include "periodic_boundaries.c"

double LENNAR_JONES_CUT_OFF_IN_SIGMA_UNIT = 2.5;
double LENNAR_JONES_LOW_DISTANCE_CUTOFF = 1e-8;

/**
 * @brief Compute the lennar jones potential associated with a particle interacting
 * with all the others, in periodic boundary condition.
 * It work regarding the space dimension.
 *
 * Apply space cutoff r_c using the fact that the lennar jones potential
 * is low range and a lot of interaction can be truncated without
 * losing accuracy.
 */
double lj_i_energy(int i,
                   const double *pos_array,
                   const double *charge_array,
                   int n_particles,
                   int space_dim,
                   double box_size,
                   double epsilon,
                   double sigma)
{
    double energy_i = 0;
    double sigma_6 = sigma * sigma * sigma * sigma * sigma * sigma;
    double sigma_12 = sigma_6 * sigma_6;

    double r_c = LENNAR_JONES_CUT_OFF_IN_SIGMA_UNIT * sigma;

    // static double VSHIFT = 0;
    // if (VSHIFT == 0)
    // {
    //     double inv_rc2 = 1. / (r_c * r_c);
    //     double inv_rc6 = inv_rc2 * inv_rc2 * inv_rc2;
    //     double inv_rc12 = inv_rc6 * inv_rc6;
    //     VSHIFT = 4.0 * epsilon * (sigma_12 * inv_rc12 - sigma_6 * inv_rc6);
    // }

    for (int k = 0; k < n_particles; k++)
    {

        // Avoid self interaction
        if (i == k)
        {
            continue;
        }

        double r2 = 0.0;

        for (int j = 0; j < space_dim; j++)
        {
            double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
            r2 += dx * dx;
        }

        if (r2 > r_c * r_c)
        {
            continue;
        }

        // Low cutoff in order to avoid computation error
        if (r2 < LENNAR_JONES_LOW_DISTANCE_CUTOFF)
        {
            r2 = LENNAR_JONES_LOW_DISTANCE_CUTOFF;
        }

        double inv_r2 = 1. / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;

        double V_Lennar_Jones = 4.0 * epsilon * (sigma_12 * inv_r12 - sigma_6 * inv_r6);

        // energy_i += V_Lennar_Jones - VSHIFT;
        energy_i += V_Lennar_Jones;
    }

    return energy_i;
}

/**
 * @brief Compute the total lennar jones potential in periodic boundary condition using "lj_i_energy".
 * See that for more.
 * It work regarding the space dimension.
 */
double lj_total_energy(const double *pos_array,
                       const double *charge_array,
                       int n_particles,
                       int space_dim,
                       double box_size,
                       double epsilon,
                       double sigma)
{
    double energy = 0.0;

    for (size_t i = 0; i < n_particles; i++)
    {
        energy += lj_i_energy(i, pos_array, charge_array, n_particles, space_dim, box_size, epsilon, sigma);
    }

    energy *= 0.5; // remove double counting from pb_compute_one_particle_energy

    return energy;
}

double lj_verlet_i_energy(int i,
                          const double *pos_array,
                          const double *charge_array,
                          const IndexesList_t *vl,
                          int n_particles,
                          int space_dim,
                          double box_size,
                          double epsilon,
                          double sigma)
{
    double energy_i = 0;

    double sigma_6 = sigma * sigma * sigma * sigma * sigma * sigma;
    double sigma_12 = sigma_6 * sigma_6;

    /**
     * Apply space cutoff using the fact that the lennar jones potential
     * is low range and a lot of interaction can be truncated without
     * losing accuracy.
     */
    double r_c = LENNAR_JONES_CUT_OFF_IN_SIGMA_UNIT * sigma;

    // static double VSHIFT = 0;
    // if (VSHIFT == 0)
    // {
    //     double inv_rc2 = 1. / (r_c * r_c);
    //     double inv_rc6 = inv_rc2 * inv_rc2 * inv_rc2;
    //     double inv_rc12 = inv_rc6 * inv_rc6;
    //     VSHIFT = 4.0 * epsilon * (sigma_12 * inv_rc12 - sigma_6 * inv_rc6);
    // }

    for (int v = 0; v < vl[i].count; v++)
    {

        size_t k = vl[i].list[v];

        // Avoid self interaction, should not be a problem
        // the verlet list does not contain self interaction
        if (i == k)
        {
            continue;
        }

        double r2 = 0.0;

        for (int j = 0; j < space_dim; j++)
        {
            double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
            r2 += dx * dx;
        }

        if (r2 > r_c * r_c)
        {
            continue;
        }

        // Low cutoff in order to avoid computation error
        if (r2 < LENNAR_JONES_LOW_DISTANCE_CUTOFF)
        {
            r2 = LENNAR_JONES_LOW_DISTANCE_CUTOFF;
        }

        double inv_r2 = 1. / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;

        double V_Lennar_Jones = 4.0 * epsilon * (sigma_12 * inv_r12 - sigma_6 * inv_r6);

        // energy_i += V_Lennar_Jones - VSHIFT;
        energy_i += V_Lennar_Jones;
    }

    return energy_i;
}

double lj_verlet_total_energy(const double *pos_array,
                              const double *charge_array,
                              const IndexesList_t *vl,
                              int n_particles,
                              int space_dim,
                              double box_size,
                              double epsilon,
                              double sigma)
{
    double energy = 0.0;

    for (size_t i = 0; i < n_particles; i++)
    {
        energy += lj_verlet_i_energy(i, pos_array, charge_array, vl, n_particles, space_dim, box_size, epsilon, sigma);
    }

    energy *= 0.5; // remove double counting from pb_compute_one_particle_energy

    return energy;
}

double lj_tail_correction_per_particle(
    double density,
    double epsilon,
    double sigma,
    double r_c)
{
    double rc3 = pow(sigma / r_c, 3);
    double rc9 = rc3 * rc3 * rc3;

    return (8.0 * PI * density * epsilon / 3.0) *
           ((1.0 / 3.0) * rc9 - rc3);
}

#endif