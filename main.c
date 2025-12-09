#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/structures.c"
#include "src/constants.c"

void save_position_state(FILE* file, float* pos_array, int n_particles, int space_dim){

    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dim; j++)
        {
            fprintf(file, "%f;", pos_array[c(i,j)]);   
        }
        fprintf(file, "\n");
    }
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
            pos_array[c(i,j)] = rand() / (RAND_MAX - 1.) * BOX_SIZE;
            vel_array[c(i,j)] = 0;
        }

        charge_array[i] = 1;
        mass_array[i] = 1;
    }
    
    clock_t begin = clock();

    /*
        EVERYTHING WE WILL WRITE 
    */
    
    clock_t end = clock();

    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;

    printf("Total time: %.0lf ms", time_spent * 1000);

    FILE* position_file = fopen("./output/position_file.csv", "w");

    save_position_state(position_file, pos_array, N, SPACE_DIM);

    free(pos_array);
    free(vel_array);
    free(mass_array);
    free(charge_array);

    fclose(position_file);

    return 0;
}
