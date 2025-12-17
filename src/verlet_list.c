#ifndef VERLET_LIST_H
#define VERLET_LIST_H

// Expected to be desnity * sphere volume
#define MAX_NEIGHBORS 400

#include <stdlib.h>
#include <stdio.h>

#include "constants.c"

struct VerletList_t
{
    int count;
    int list[MAX_NEIGHBORS]; // List of indexes
} typedef VerletList_t;

void verlet_pb_build_list(const double *pos_array,
                          double *old_pos_array,
                          VerletList_t *vl,
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
                if (vl[i].count < MAX_NEIGHBORS)
                {
                    vl[i].list[vl[i].count++] = k;
                }
                else
                {
                    fprintf(stderr, "ERROR: MAX_NEIGHBORS exceeded for particle %zu\n", i);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    memcpy(old_pos_array, pos_array, total_array_size);
}

int verlet_pb_needs_rebuild(const double *pos_array,
                            const double *old_pos_array,
                            int n_particles,
                            int space_dim,
                            double box_size,
                            double skin)
{
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

int get_max_verlet_count(VerletList_t *vl,
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