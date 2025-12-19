#ifndef CHECKPOINTS_HANDLER_C
#define CHECK_POINT_HANDLER_C

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Function to save the current system and rng state in a bin file.
 * Used to start a simulation from a fixed point in time instead of having to simulate all over again.
 * Usefull for study more precisely certain range of parameters.
 */
void save_checkpoint_binary(const char *filename,
                            const double *pos_array,
                            int n_particles,
                            int space_dim,
                            double energy)
{
    FILE *f = fopen(filename, "wb");
    if (!f)
    {
        perror("Checkpoint save failed");
        exit(EXIT_FAILURE);
    }

    fwrite(&n_particles, sizeof(int), 1, f);
    fwrite(&space_dim, sizeof(int), 1, f);

    fwrite(&energy, sizeof(double), 1, f);

    fwrite(pos_array, sizeof(double), n_particles * space_dim, f);

    // Save rng state
    unsigned short rng_state[3];
    erand48(rng_state);
    fwrite(rng_state, sizeof(unsigned short), 3, f);

    fclose(f);
}

/**
 * @brief Function to load a system and rng state from a bin file.
 * Used to start a simulation from a fixed point in time instead of having to simulate all over again.
 * Usefull for study more precisely certain range of parameters.
 */
void load_checkpoint_binary(const char *filename,
                            double *pos_array,
                            int n_particles,
                            int space_dim,
                            double *energy)
{
    FILE *f = fopen(filename, "rb");
    if (!f)
    {
        perror("Checkpoint load failed");
        exit(EXIT_FAILURE);
    }

    int n_p, s_d;
    fread(&n_p, sizeof(int), 1, f);
    fread(&s_d, sizeof(int), 1, f);

    if (n_p != n_particles || s_d != space_dim)
    {
        fprintf(stderr, "ERROR: Incompatible checkpoint. N particle or space dimension missmatch\n");
        exit(EXIT_FAILURE);
    }

    fread(energy, sizeof(double), 1, f);

    fread(pos_array, sizeof(double),
          n_particles * space_dim, f);

    // Load rng state
    unsigned short rng_state[3];
    fread(rng_state, sizeof(unsigned short), 3, f);
    srand48(rng_state[0]);

    fclose(f);
}

#endif