#ifndef CHECK_POINTS_HANDLER_C
#define CHECK_POINTS_HANDLER_C

#include <stdio.h>
#include <stdlib.h>

#include "logger.c"

static const char CHECKPOINT_ERROR_MESSAGE_INFO_MISSMATCH[] = "Incompatible checkpoint. N particle or space dimension missmatch";
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
        LOG_FATAL("Checkpoint save failed");
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
        LOG_FATAL("Checkpoint load failed");
    }

    int n_p, s_d;
    size_t chunk_return;

    chunk_return = fread(&n_p, sizeof(int), 1, f);
    chunk_return = fread(&s_d, sizeof(int), 1, f);

    if (n_p != n_particles || s_d != space_dim)
    {
        LOG_FATAL("%s", CHECKPOINT_ERROR_MESSAGE_INFO_MISSMATCH);
    }

    chunk_return = fread(energy, sizeof(double), 1, f);
    chunk_return = fread(pos_array, sizeof(double),
                         n_particles * space_dim, f);

    // Load rng state
    unsigned short rng_state[3];
    chunk_return = fread(rng_state, sizeof(unsigned short), 3, f);
    srand48(rng_state[0]);

    if (chunk_return == 1)
        printf("here to remove a warming, just ignore");

    fclose(f);
}

#endif
