#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/constants.c"
#include "src/structures.c"

void save_position_state(FILE *file, float *pos_array, int n_particles, int space_dim)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dim; j++)
        {
            fprintf(file, "%f;", pos_array[c(i, j)]);
        }
        fprintf(file, "\n");
    }
}

float compute_energy(float *pos_array, float *charge_array, int n_particles, int space_dim)
{
    float energy = 0;

    for (size_t i = 0; i < n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < n_particles; j++)
        {

            float distance_square = 0;

            for (size_t k = 0; k < space_dim; k++)
            {
                distance_square += (pos_array[c(i, k)] - pos_array[c(j, k)]) * (pos_array[c(i, k)] - pos_array[c(j, k)]);
            }

            if (distance_square == 0)
            {
                printf("WARNING: Overlapping particle found (i=%ld,j=%ld)\n", i, j);
                continue;
            }

            energy += (charge_array[i] * charge_array[j]) / sqrtf(distance_square);
        }
    }

    return energy;
}

int main(int argc, char const *argv[])
{
    // Set seed for reproducibility
    srand(SEED);

    int total_vel_pos_array_size = N * SPACE_DIM;

    // Positions and velocities are store in an array of size (N * SPACE_DIM) where
    // the j component of the i particle is stored at index (SPACE_DIM * i + j).
    // In order to retrieve the correct index for component j of particle i use the macro "c(i,j)"

    float *pos_array = (float *)malloc(total_vel_pos_array_size * sizeof(float));
    if (pos_array == NULL)
        exit(EXIT_FAILURE);

    float *vel_array = (float *)malloc(total_vel_pos_array_size * sizeof(float));
    if (vel_array == NULL)
        exit(EXIT_FAILURE);

    float *charge_array = (float *)malloc(N * sizeof(float));
    if (charge_array == NULL)
        exit(EXIT_FAILURE);

    float *mass_array = (float *)malloc(N * sizeof(float));
    if (mass_array == NULL)
        exit(EXIT_FAILURE);

        
    // Init arrays
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < SPACE_DIM; j++)
        {
            // Uniform position distribution inside the square box
            pos_array[c(i, j)] = rand() / (RAND_MAX - 1.) * BOX_SIZE;
            vel_array[c(i, j)] = 0;
        }

        charge_array[i] = 1;
        mass_array[i] = 1;
    }

    // Evaluate the execution time
    clock_t begin = clock();

    /*
        EVERYTHING WE WILL WRITE
    */

    clock_t end = clock();

    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;

    printf("Total time: %.0lf ms\n", time_spent * 1000);

    FILE *position_file = fopen("./output/position_file.csv", "w");
    FILE *results = fopen("./output/results.csv", "w");

    save_position_state(position_file, pos_array, N, SPACE_DIM);

    printf("Coulomb energy: %f\n", compute_energy(pos_array, charge_array, N, SPACE_DIM));
    
    free(pos_array);
    free(vel_array);
    free(mass_array);
    free(charge_array);

    fclose(position_file);

    return 0;
}
