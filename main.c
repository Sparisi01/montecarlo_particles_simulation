#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "src/constants.c"
#include "src/structures.c"

// Generate velocities from boltzman distribution using box muller
double getRNDVelocity(double temperature)
{
    double x = drand48();
    double y = drand48();

    return sqrt(temperature) * sqrt(-2 * log(1 - x)) * cos(2 * PI * y);
}

void save_particle_state(FILE *file, double *pos_array, double *charge_array, int n_particles, int space_dim)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dim; j++)
        {
            fprintf(file, "%f;", pos_array[c(i, j)]);
        }
        fprintf(file, "%f\n", charge_array[i]);
    }
}

static inline double pair_potential(double qi, double qj, double r2)
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
    //double V_Lennar_Jones = 4.0 * EPSILON * (sigma_12 * inv_r12 - sigma_6 * inv_r6);
    double V_quadratic = inv_r4;
    
    return V_coulomb + V_quadratic;
}

/* COMPUTE_ONE_PARTICLE_ENERGY
 *
 * Computes the interaction energy of a single particle `i` with all other particles
 * in the system. The total energy is obtained by summing the pairwise potential
 * between particle `i` and each particle `j`..
 *
 * Parameters:
 *   i             - Index of the particle for which the energy is computed.
 *   pos_array     - Array containing the positions of all particles. Coordinates
 *                   are accessed using the macro c(i, k).
 *   charge_array  - Array containing the charges of all particles, used for Coulomb potential.
 *   n_particles   - Total number of particles in the system.
 *   space_dim     - Dimensionality of the space (e.g., 2 for 2D, 3 for 3D).
 *
 * Behavior:
 *   - Skips self-interaction (when j == i).
 *   - Computes the squared distance between particles i and j by summing the
 *     squared coordinate differences.
 *   - Issues a warning and skips the pair if the distance is below a small
 *     threshold, indicating overlapping particles.
 *   - Accumulates the energy contribution from the pair potential function.
 *
 * Returns:
 *   The total interaction energy of particle `i` as a double.
 */

double compute_one_particle_energy(int i, double *pos_array, double *charge_array, int n_particles, int space_dim)
{
    double E = 0;

    // #pragma omp parallel for reduction(+:energy)
    for (int j = 0; j < n_particles; j++)
    {   
        // Avoid self interaction
        if (i == j)
            continue;

        double r2 = 0;

        for (int k = 0; k < space_dim; k++)
        {
            double dx = pos_array[c(i, k)] - pos_array[c(j, k)];
            r2 += dx * dx;
        }

        E += pair_potential(charge_array[i], charge_array[j], r2);
    }

    return E;
}

// Convert the energy array to the total energy
double array_to_total_energy(double *E_array, int n_particles)
{
    double E = 0;

    for (int i = 0; i < n_particles; i++)
    {
        E += E_array[i];
    }

    return E / 2;
}

double compute_average_pair_distance(double *pos_array, double *charge_array, int n_particles, int space_dim)
{
    double sum_pair_distance = 0;
    FILE *file_distances = fopen("./output/file_distances.csv", "w");

    for (int i = 0; i < n_particles - 1; i++)
    {

        for (int j = i + 1; j < n_particles; j++)
        {
            if (i == j)
                continue;

            double distance_square = 0;

            for (int k = 0; k < space_dim; k++)
            {
                distance_square += (pos_array[c(i, k)] - pos_array[c(j, k)]) * (pos_array[c(i, k)] - pos_array[c(j, k)]);
            }

            if (distance_square < 1e-8)
            {
                printf("\r\033[2K WARNING: Overlapping particle found (i=%d,j=%d)\n", i, j);
                continue;
            }
            if ((charge_array[i] == 4 && charge_array[j] == -1) || (charge_array[j] == 4 && charge_array[i] == -1))
                fprintf(file_distances, "%lf\n", sqrt(distance_square));
            sum_pair_distance += sqrt(distance_square);
        }
    }

    return (sum_pair_distance * 2) / (N * (N - 1.));
}

void metropolis_step(double *energy_array, double *pos_array, double *charge_array, double delta, double temperature, int n_particles, int space_dim, long *accepted_counter)
{
    // Alocate an array of dj on the stack,
    // this is used to keep track of the different direction and make the code independent on space dimension
    double dj_array[space_dim];

    // Update one particle at the time. Better to check box boundaries
    for (int i = 0; i < n_particles; i++)
    {

        // NOTE: right now I am computing old and new energy, each is O(N)
        double old_i_energy = compute_one_particle_energy(i, pos_array, charge_array, n_particles, space_dim);

        for (int j = 0; j < space_dim; j++)
        {
            // Random step in j direction between -delta and + delta
            double dj = ((drand48() * 2) - 1) * delta;

            // Refuse step if particle i in direction j out of the box
            if (pos_array[c(i, j)] + dj > BOX_SIZE || pos_array[c(i, j)] + dj < 0)
            {
                dj = 0;
            }

            pos_array[c(i, j)] += dj;
            dj_array[j] = dj;
        }

        /** NOTE: Since I am using compute_one_particle_energy I can not use the fact that the interaction between i and k
         * is the same as k and i in order to cut the computation in half. Nevertheless it is usefull because in this
         * way I can check out of the boundaries conditions for one particle at the time,
         * instead of having to discard a full system step just because one particle
         * felt off.
         */
        double new_i_energy = compute_one_particle_energy(i, pos_array, charge_array, n_particles, space_dim);

        //double dE = new_i_energy - energy_array[i];
        double dE = new_i_energy - old_i_energy;

        // METROPOLIS ACCEPTANCE AND UPDATE energy_array
        int accepted = 0;
        double alpha = fmin(1, exp( -dE / temperature));
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
            // Update i_particle energy if step accepted
            energy_array[i] = new_i_energy;
            (*accepted_counter)++;
        }
    }
}

void init_system(double *pos_array, double *charge_array, double *mass_array)
{
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < SPACE_DIM; j++)
        {
            // Uniform position distribution inside the square box
            pos_array[c(i, j)] = drand48() * BOX_SIZE;
            // vel_array[c(i, j)] = getRNDVelocity(1);
        }

        switch (i % 5)
        {
        case 0:
            charge_array[i] = 4;
            break;

        default:
            charge_array[i] = -1;
            break;
        }

        mass_array[i] = 1;
    }
}

/* PRINT_PROGRESS
 *
 * Prints a terminal progress bar with percentage and estimated time remaining (ETA).
 *
 * Parameters:
 *   current     – current iteration (0 … total)
 *   total       – total number of iterations
 *   start_time  – the time (time_t) when the process started; used to compute ETA
 *
 * Behavior:
 *   - The function prints a 50-character progress bar using full block characters (█)
 *     and a partial block cursor (▌).
 *   - It overwrites the same terminal line using '\r'.
 *   - ETA is computed from elapsed time and current progress and displayed as HH:MM:SS.
 *
 * Usage:
 *   Call this function periodically inside a loop to update the progress bar without
 *   significantly impacting performance.
 */
void print_progress(size_t current, size_t total, clock_t start_time)
{
    const int barWidth = 50;
    double progress = (double)current / total;
    int filled = progress * barWidth;

    // Time elapsed
    clock_t now = clock();
    double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;

    // ETA stimation
    double eta = (progress > 0.0) ? elapsed * (1.0 - progress) / progress : 0.0;

    // Conversion of ETA hh:mm:ss
    int eta_h = (int)eta / 3600;
    int eta_m = ((int)eta % 3600) / 60;
    int eta_s = (int)eta % 60;

    printf("\r[");

    for (int i = 0; i < barWidth; i++)
    {
        if (i < filled)
            printf("█");
        else if (i == filled)
            printf("▌");
        else
            printf(" ");
    }

    printf("] %5.1f%%  ETA: %02d:%02d:%02d (hh:mm:ss)",
           progress * 100, eta_h, eta_m, eta_s);

    fflush(stdout);
}

int main(int argc, char const *argv[])
{
    // omp_set_num_threads(8);
    //  Set seed for reproducibility
    srand48(SEED);

    int total_vel_pos_array_size = N * SPACE_DIM;


    //-------------------------------------------
    // Positions and velocities are store in an array of size (N * SPACE_DIM) where
    // the j component of the i particle is stored at index (SPACE_DIM * i + j).
    // In order to retrieve the correct index for component j of particle i use the macro "c(i,j)"

    double *pos_array = (double *)malloc(total_vel_pos_array_size * sizeof(double));
    if (pos_array == NULL)
        exit(EXIT_FAILURE);

    double *vel_array = (double *)malloc(total_vel_pos_array_size * sizeof(double));
    if (vel_array == NULL)
        exit(EXIT_FAILURE);

    double *charge_array = (double *)malloc(N * sizeof(double));
    if (charge_array == NULL)
        exit(EXIT_FAILURE);

    double *mass_array = (double *)malloc(N * sizeof(double));
    if (mass_array == NULL)
        exit(EXIT_FAILURE);

    double *energy_array = (double *)malloc(N * sizeof(double));
    if (energy_array == NULL)
        exit(EXIT_FAILURE);

    FILE *start_position_file = fopen("./output/start_position_file.csv", "w");
    if (start_position_file == NULL)
        exit(EXIT_FAILURE);
    
    FILE *end_position_file = fopen("./output/end_position_file.csv", "w");
    if (end_position_file == NULL)
        exit(EXIT_FAILURE);

    FILE *energy_file = fopen("./output/energy.csv", "w");
    if (energy_file == NULL)
        exit(EXIT_FAILURE);

    //-------------------------------------------
    
    // Init system
    init_system(pos_array, charge_array, mass_array);

    // Init energy array
    for (size_t i = 0; i < N; i++)
    {
        energy_array[i] = compute_one_particle_energy(i, pos_array, charge_array, N, SPACE_DIM);
    }

    save_particle_state(start_position_file, pos_array, charge_array, N, SPACE_DIM);
    printf("Start energy: %f\n", array_to_total_energy(energy_array, N));

    //----------------------------------
    // START Evaluate the execution time
    clock_t begin = clock();
    //----------------------------------

    long accepted_steps = 0;

    for (size_t i = 0; i < N_METROPOLIS_STEPS; i++)
    {
        metropolis_step(energy_array, pos_array, charge_array, STEP_SIZE, TEMPERATURE, N, SPACE_DIM, &accepted_steps);
        fprintf(energy_file, "%lf\n", array_to_total_energy(energy_array, N));

        // Progress Bar
        if (i % PRINT_INTERVAL == 0)
            print_progress(i, N_METROPOLIS_STEPS, begin);
    }

    // Clear terminal
    printf("\r\033[2K");

    //----------------------------------
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time: %.0lf ms\n", time_spent * 1000);
    // END total time evaliation
    //----------------------------------

    printf("Accepted step: %lf\n", accepted_steps / ((double)N_METROPOLIS_STEPS * N));
    printf("Average pair distance: %lf\n", compute_average_pair_distance(pos_array, charge_array, N, SPACE_DIM));
    printf("End energy: %f\n", array_to_total_energy(energy_array, N));

    save_particle_state(end_position_file, pos_array, charge_array, N, SPACE_DIM);

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
