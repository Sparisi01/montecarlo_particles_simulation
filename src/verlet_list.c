#ifndef VERLET_LIST_C
#define VERLET_LIST_C

// Expected to be density * sphere volume
#define VERLET_MAX_NEIGHBORS 1024

#include <stdlib.h>
#include <stdio.h>

#include "constants.c"
#include "periodic_boundaries.c"
#include "logger.c"

static const char VERLET_ERROR_MAX_NEIGHT_EXCEED[] = "Max neightbors exceeded.\n"
                                                     "Increase VERLET_MAX_NEIGHBORS.";

struct IndexesList_t
{
    int count;
    int list[VERLET_MAX_NEIGHBORS]; // List of indexes
} typedef IndexesList_t;

void verlet_build_list(const double *pos_array,
                       double *old_pos_array,
                       IndexesList_t *vl,
                       int n_particles,
                       int space_dim,
                       double box_size,
                       double r_cut,
                       double skin)
{
    const double r_verlet2 = (r_cut + skin) * (r_cut + skin);
    const size_t total_array_size = (size_t)n_particles * space_dim * sizeof(double);

    for (size_t i = 0; i < n_particles; i++)
    {
        vl[i].count = 0;

        for (size_t k = 0; k < n_particles; k++)
        {

            if (k == i)
                continue;

            double r2 = 0.0;
            for (size_t j = 0; j < space_dim; j++)
            {
                double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
                r2 += dx * dx;
            }

            if (r2 < r_verlet2)
            {
                if (vl[i].count < VERLET_MAX_NEIGHBORS)
                {
                    vl[i].list[vl[i].count++] = k;
                }
                else
                {
                    LOG_FATAL("%s", VERLET_ERROR_MAX_NEIGHT_EXCEED);
                }
            }
        }
    }

    memcpy(old_pos_array, pos_array, total_array_size);
}

int verlet_check_needs_rebuild(const double *pos_array,
                               const double *old_pos_array,
                               int n_particles,
                               int space_dim,
                               double box_size,
                               double skin)
{
    /** NOTE Multiplied by 0.5 to stay safe in the case that two particles moved toward each other
     *   by (skin * 0.5), resulting in a new displacement distance decreased by (2*skin*0.5) = skin
     */
    double limit2 = (skin * 0.5) * (skin * 0.5);

    for (size_t i = 0; i < n_particles; i++)
    {

        double r2 = 0.0;
        for (size_t j = 0; j < space_dim; j++)
        {
            double dx = pb_minimum_image(pos_array[c(i, j)] - old_pos_array[c(i, j)], box_size);
            r2 += dx * dx;
        }

        if (r2 > limit2)
            return 1;
    }
    return 0;
}

/**
 * @brief Get max number of neightbours in the verlet list
 */
int verlet_get_max_neightbours(IndexesList_t *vl,
                               int n_particles)
{
    int max_count = 0;

    for (size_t i = 0; i < n_particles; i++)
    {
        if (vl[i].count > max_count)
            max_count = vl[i].count;
    }
    return max_count;
}
#endif