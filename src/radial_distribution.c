#ifndef RADIAL_DISTRIBUTION_C
#define RADIAL_DISTRIBUTION_C

#include <stdlib.h>

/**
 * @brief Compute the radial distribution storing the counting in a bin array "bins_array" of bin size "bin_interval"
 */
void radial_distribution(const double *pos_array,
                         int n_particles,
                         int space_dim,
                         double box_size,
                         double *bins_array,
                         int N_bins,
                         double bin_interval)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t k = i + 1; k < n_particles; k++)
        {
            double r2 = 0;
            for (int j = 0; j < space_dim; j++)
            {
                double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
                r2 += dx * dx;
            }

            int n_bin = (int)(sqrt(r2) / bin_interval);

            if (n_bin < N_bins)
                bins_array[n_bin]++;
        }
    }
}

void radial_distribution_diff_charges(const double *pos_array,
                                      const double *charge_array,
                                      int n_particles,
                                      int space_dim,
                                      double box_size,
                                      double *bins_array,
                                      int N_bins,
                                      double bin_interval)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t k = i + 1; k < n_particles; k++)
        {
            // Keep track only of particle with different charges
            if (charge_array[i] == charge_array[k])
                continue;

            double r2 = 0;
            for (int j = 0; j < space_dim; j++)
            {
                double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
                r2 += dx * dx;
            }

            int n_bin = (int)(sqrt(r2) / bin_interval);

            if (n_bin < N_bins)
                bins_array[n_bin]++;
        }
    }
}

/**
 * @brief Compute the radial distribution storing the counting in a bin array
 */
void radial_distribution_equal_charges(const double *pos_array,
                                       const double *charge_array,
                                       int n_particles,
                                       int space_dim,
                                       double box_size,
                                       double *bins_array,
                                       int N_bins,
                                       double bin_interval)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t k = i + 1; k < n_particles; k++)
        {
            // Keep track only of particle with equal charges
            if (charge_array[i] != charge_array[k])
                continue;

            double r2 = 0;
            for (int j = 0; j < space_dim; j++)
            {
                double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
                r2 += dx * dx;
            }

            int n_bin = (int)(sqrt(r2) / bin_interval);

            if (n_bin < N_bins)
                bins_array[n_bin]++;
        }
    }
}

// Compute all three possible radial distribution at once using only one for loop
void radial_distribution_all(const double *pos_array,
                             const double *charge_array,
                             int n_particles,
                             int space_dim,
                             double box_size,
                             double *bins_array,
                             double *bins_array_equal,
                             double *bins_array_differ,
                             int N_bins,
                             double bin_interval)
{

    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t k = i + 1; k < n_particles; k++)
        {

            double r2 = 0;
            for (int j = 0; j < space_dim; j++)
            {
                double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
                r2 += dx * dx;
            }

            int n_bin = (int)(sqrt(r2) / bin_interval);

            if (n_bin < N_bins)
            {

                if (charge_array[i] == charge_array[k])
                    bins_array_equal[n_bin]++;

                if (charge_array[i] != charge_array[k])
                    bins_array_differ[n_bin]++;

                bins_array[n_bin]++;
            }
        }
    }
}

#endif