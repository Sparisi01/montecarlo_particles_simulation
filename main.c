#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/constants.c"
#include "src/structures.c"

// Generate velocities from boltzman distribution
// using box muller
float getRNDVelocity(float temperature)
{
    double x = rand() / (RAND_MAX + 1.);
    double y = rand() / (RAND_MAX + 1.);

    return sqrtf(temperature) * sqrt(-2 * log(1 - x)) * cos(2 * PI * y);
}

void save_position_state(FILE *file, float *pos_array, int n_particles, int space_dim)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dim; j++)
        {
            fprintf(file, "%f", pos_array[c(i, j)]);

            if(j != space_dim - 1)
                fprintf(file, ";");
        }
        fprintf(file, "\n");
    }
}

float compute_energy(float *pos_array, float *charge_array, int n_particles, int space_dim)
{
    float energy = 0;

    for (int i = 0; i < n_particles - 1; i++)
    {
        for (int j = i + 1; j < n_particles; j++)
        {
            float distance_square = 0;

            for (int k = 0; k < space_dim; k++)
            {
                distance_square += (pos_array[c(i, k)] - pos_array[c(j, k)]) * (pos_array[c(i, k)] - pos_array[c(j, k)]);
            }

            if (distance_square == 0)
            {
                printf("WARNING: Overlapping particle found (i=%d,j=%d)\n", i, j);
                continue;
            }

            energy += (charge_array[i] * charge_array[j]) / sqrtf(distance_square);
        }
    }
    return energy;
}

float compute_one_particle_energy(int i, float *pos_array, float *charge_array, int n_particles, int space_dim)
{
    float energy = 0;

    for (int j = 0; j < n_particles; j++)
    {
        if (i == j) continue;

        float distance_square = 0;

        for (int k = 0; k < space_dim; k++)
        {
            distance_square += (pos_array[c(i, k)] - pos_array[c(j, k)]) * (pos_array[c(i, k)] - pos_array[c(j, k)]);
        }

        if (distance_square == 0)
        {
            printf("WARNING: Overlapping particle found (i=%d,j=%d)\n", i, j);
            continue;
        }

        energy += (charge_array[i] * charge_array[j]) / sqrtf(distance_square); // Coulomb
        energy += 1/pow(distance_square,4); // Hard core
    }

    return energy;
}

// Convert the energy array to the total energy
float array_to_total_energy(float *energy_array, int n_particles)
{
    float energy = 0;

    for(int i = 0; i<n_particles; i++){
        energy += energy_array[i];
    }

    return energy/2;
}

float metropolis_step(float *energy_array, float *pos_array, float *charge_array, float delta, float temperature, int n_particles, int space_dim)
{
    // Alocate an array of dj on the stack,
    // this is used to keep track of the different direction
    float dj_array[space_dim];

    // Update one particle at the time. Better to check box boundaries
    for (int i = 0; i < n_particles; i++)
    {
        for (int j = 0; j < space_dim; j++)
        {
            // Random step in j direction between -delta and + delta
            float dj = (((rand() / (RAND_MAX + 1.))*2)-1)*delta;
            
            // Refuse step if out of the box
            if (pos_array[c(i, j)] + dj > BOX_SIZE || pos_array[c(i, j)] + dj < 0)
            {   
                dj = 0;
            }

            pos_array[c(i, j)] += dj;
            dj_array[j] = dj;
        }

        float new_i_energy = compute_one_particle_energy(i, pos_array, charge_array, n_particles, space_dim);

        // METROPOLIS ACCEPTANCE AND UPDATE energy_array

        int accepted = 0;
        float alpha = fmin(1, exp((energy_array[i] - new_i_energy) / temperature));
        accepted = rand() / (RAND_MAX + 1.) <= alpha;

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
            // Update particle i energy if step accepted
            energy_array[i] = new_i_energy;
        }
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

    float *energy_array = (float *)malloc(N * sizeof(float));
    if (energy_array == NULL)
        exit(EXIT_FAILURE);
    
    // Init arrays
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < SPACE_DIM; j++)
        {
            // Uniform position distribution inside the square box
            pos_array[c(i, j)] = rand() / (RAND_MAX + 1.) * BOX_SIZE;
            // vel_array[c(i, j)] = getRNDVelocity(1);
        }

        charge_array[i] = rand() / (RAND_MAX + 1.) > 0.5;
        mass_array[i] = 1;
    }

    // Init energy array
    for (size_t i = 0; i < N; i++)
    {
        energy_array[i] = compute_one_particle_energy(i, pos_array, charge_array, N, SPACE_DIM);
    }

    // Save starting particle position
    FILE *start_position_file = fopen("./output/start_position_file.csv", "w");
    save_position_state(start_position_file, pos_array, N, SPACE_DIM);

    FILE *energy_file = fopen("./output/energy.csv", "w");

    printf("Start energy: %f\n", array_to_total_energy(energy_array, N));
    
    // START Evaluate the execution time
    clock_t begin = clock();

    

    for (size_t i = 0; i < N_METROPOLIS_STEPS; i++)
    {
        metropolis_step(energy_array, pos_array, charge_array, 1e-1, 1, N, SPACE_DIM);
        fprintf(energy_file, "%lf\n", array_to_total_energy(energy_array, N));
    }

   

    clock_t end = clock();

    // END total time evaliation

    printf("End energy: %f\n", array_to_total_energy(energy_array, N));

    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time: %.0lf ms\n", time_spent * 1000);

    // Save ending particle position
    FILE *end_position_file = fopen("./output/end_position_file.csv", "w");
    save_position_state(end_position_file, pos_array, N, SPACE_DIM);
    

    free(pos_array);
    free(vel_array);
    free(mass_array);
    free(charge_array);
    free(energy_array);

    fclose(start_position_file);
    fclose(end_position_file);
    fclose(energy_file);

    return 0;
}
