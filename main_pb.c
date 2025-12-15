#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>

#include "src/constants.c"
#include "src/periodic_boundaries.c"
#include "src/ewald.c"
#include "src/progress_bar.c"

double array_mean(double *arr, int n)
{
    double mean = 0;
    for (size_t i = 0; i < n; i++)
    {
        mean += arr[i];
    }
    return mean / n;
}

double array_mean2(double *arr, int n)
{
    double mean2 = 0;
    for (size_t i = 0; i < n; i++)
    {
        mean2 += arr[i] * arr[i];
    }
    return mean2 / n;
}

double array_var(double *arr, int n)
{
    return array_mean2(arr, n) - array_mean(arr, n) * array_mean(arr, n);
}

void radial_distribution(const double *pos_array, int n_particles, int space_dim, double box_size, double *bins_array, int N_bins, double bin_interval)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t k = 0; k < n_particles; k++)
        {
            if (k == i)
                continue;

            double r2 = 0;
            for (int j = 0; j < space_dim; j++)
            {
                double dx = pb_minimum_image(pos_array[c(i, j)] - pos_array[c(k, j)], box_size);
                r2 += dx * dx;
            }

            int n_bin = (int)floor(sqrt(r2) / bin_interval);
            bins_array[n_bin]++;
        }
    }
}

void save_particle_state(FILE *file, double *pos_array, double *charge_array, int n_particles, int space_dim)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dim; j++)
        {
            double x = pos_array[c(i, j)];
            fprintf(file, "%f;", pos_array[c(i, j)]);
        }
        fprintf(file, "%f\n", charge_array[i]);
    }
}

double compute_i_lennar_jones_potential(int i, const double *pos_array, const double *charge_array, int n_particles, int space_dim, double box_size)
{
    double energy_i = 0;

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

        double r_c = 2.5 * SIGMA;

        // Apply cutoff
        if (r2 > r_c * r_c)
        {
            continue;
        }

        // Low cutoff in order to avoid computation error
        if (r2 < 1e-6)
        {
            r2 = 1e-6;
        }

        double inv_r2 = 1. / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;

        double sigma_6 = SIGMA * SIGMA * SIGMA * SIGMA * SIGMA * SIGMA;
        double sigma_12 = sigma_6 * sigma_6;

        double VSHIFT = (4 * EPSILON * (pow(SIGMA / (box_size / 2), 12) - pow(SIGMA / (box_size / 2), 6)));
        double V_Lennar_Jones = 4.0 * EPSILON * (sigma_12 * inv_r12 - sigma_6 * inv_r6) - VSHIFT;
        energy_i += V_Lennar_Jones;
    }

    return energy_i;
}

double compute_lennar_jones_energy(const double *pos_array, const double *charge_array, int n_particles, int space_dim, double box_size)
{
    double energy = 0.0;

    for (size_t i = 0; i < n_particles; i++)
    {
        energy += compute_i_lennar_jones_potential(i, pos_array, charge_array, n_particles, space_dim, box_size);
    }
    // remove double counting from pb_compute_one_particle_energy
    energy *= 0.5;

    return energy;
}

// Given the system of N particle in D dimensional space compute the total interaction energy
double pb_compute_total_energy(const double *pos_array, const double *charge_array, int n_particles, int space_dim, double box_size)
{
    double total_energy = 0;
    total_energy += compute_lennar_jones_energy(pos_array, charge_array, n_particles, space_dim, box_size);
    // total_energy += ewd_total_coulomb_energy(pos_array, charge_array, n_particles, box_size);

    return total_energy;
}

// Initialize positions in a cubic lattice
void init_system_lattice(double *pos_array, double *charge_array, double *mass_array, int n_particles, double box_size, int lattice_type, int n_cell_per_row)
{
    struct Vec3
    {
        double x;
        double y;
        double z;
    } typedef Vec3;

    // lattice_type 1 CC, 2 BCC, 4 FCC

    Vec3 lattice_positions[lattice_type];

    switch (lattice_type)
    {
    case 1: // CC
        lattice_positions[0] = (Vec3){0, 0, 0};
        break;
    case 2: // BCC
        lattice_positions[0] = (Vec3){0, 0, 0};
        lattice_positions[1] = (Vec3){0.5, 0.5, 0.5};
        break;
    case 4: // FCC
        lattice_positions[0] = (Vec3){0, 0, 0};
        lattice_positions[1] = (Vec3){0, 0.5, 0.5};
        lattice_positions[2] = (Vec3){0.5, 0.5, 0};
        lattice_positions[3] = (Vec3){0.5, 0, 0.5};
        break;
    default:
        break;
    }

    double single_cell_length = (double)box_size / n_cell_per_row;
    int n_particle_placed = 0;

    for (size_t n_cell_x = 0; n_cell_x < n_cell_per_row; n_cell_x++)
    {
        for (size_t n_cell_y = 0; n_cell_y < n_cell_per_row; n_cell_y++)
        {
            for (size_t n_cell_z = 0; n_cell_z < n_cell_per_row; n_cell_z++)
            {
                for (size_t w = 0; w < lattice_type; w++)
                {
                    pos_array[c(n_particle_placed, 0)] = (n_cell_x + lattice_positions[w].x) * single_cell_length;
                    pos_array[c(n_particle_placed, 1)] = (n_cell_y + lattice_positions[w].y) * single_cell_length;
                    pos_array[c(n_particle_placed, 2)] = (n_cell_z + lattice_positions[w].z) * single_cell_length;

                    charge_array[n_particle_placed] = (n_particle_placed % 2 == 0) ? 1 : -1;
                    mass_array[n_particle_placed] = 1;

                    n_particle_placed++;
                }
            }
        }
    }
}

void init_system(double *pos_array, double *charge_array, double *mass_array, int n_particles, int space_dimension, double box_size)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dimension; j++)
        {
            // Uniform position distribution inside the square box
            pos_array[c(i, j)] = drand48() * box_size;
        }

        switch (i % 2)
        {
        case 0:
            charge_array[i] = 1;
            break;

        default:
            charge_array[i] = -1;
            break;
        }

        mass_array[i] = 1;
    }
}

// This function perform an update of all the system at once, hard to get good acceptance rate for stable configuration.
double pb_metropolis_step_full_system(double old_energy, double *pos_array, const double *charge_array, double delta, double temperature, int n_particles, int space_dim, long *accepted_counter, double box_size)
{
    // Save old position configuration
    // NOTE - This memory gets never free
    static double *old_pos_array = NULL;
    size_t total_array_size = sizeof(double) * n_particles * space_dim;

    if (old_pos_array == NULL)
    {
        old_pos_array = (double *)malloc(total_array_size);
        if (old_pos_array == NULL)
        {
            fprintf(stderr, "Memory allocation failed for old_pos_array\n");
            return NAN;
        }
    }

    // Copy positions from pos_array to old_pos_array
    memcpy(old_pos_array, pos_array, total_array_size);

    // Update system: try to move each particle
    for (int i = 0; i < n_particles; i++)
    {
        for (int j = 0; j < space_dim; j++)
        {
            // Random step in j direction between -delta and + delta
            double dj = ((drand48() * 2) - 1) * delta;

            // NOTE: In periodic boundaries conditions there is no need to check for boundaries
            pos_array[c(i, j)] += dj;
        }
    }

    // Compute the new energy
    double new_energy = pb_compute_total_energy(pos_array, charge_array, n_particles, space_dim, box_size);
    double dE = new_energy - old_energy;

    // Metropolis acceptance criterion
    int accepted = 0;
    double alpha = fmin(1, exp(-dE / temperature));
    accepted = drand48() <= alpha;

    // If the step is not accepted cancel the position update for the particle i
    if (!accepted)
    {
        // Retrieve old positions
        memcpy(pos_array, old_pos_array, total_array_size);

        return old_energy;
    }
    else
    {
        (*accepted_counter)++;

        // Bring the particle back to the first cell
        for (size_t i = 0; i < n_particles; i++)
        {
            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] = pb_wrap_position(pos_array[c(i, j)], box_size);
            }
        }

        return new_energy;
    }
}

// This function perform an update of all the system on particle at the time, easier to get good acceptance rate but
// it has highter temporal correlation
double pb_metropolis_step_one_particle(double energy, double *pos_array, const double *charge_array, double delta, double temperature, int n_particles, int space_dim, long *accepted_counter, double box_size)
{
    // Array of indexes
    static int *perm = NULL;
    if (perm == NULL)
    {
        perm = (int *)malloc(n_particles * sizeof(int));
        for (int i = 0; i < n_particles; i++)
        {
            perm[i] = i;
        }
    }

    // Generate a random permutation of indexes
    for (int i = n_particles - 1; i > 0; i--)
    {
        int j = (int)(drand48() * (i + 1));
        int tmp = perm[i];
        perm[i] = perm[j];
        perm[j] = tmp;
    }

    // Keep track of space shifts in memory in order to restore the previous position if the step is not accepted
    double steps_save[space_dim];

    // Update one particle at the time
    // NOTE: right now I am update all the particle in order, could be better to update i random particles with repetition
    // in order to remove bias. This option has to be studied
    for (size_t p = 0; p < n_particles; p++)
    {

        int i = perm[p];
        // Calculate old potential
        double old_energy = compute_i_lennar_jones_potential(i, pos_array, charge_array, n_particles, space_dim, box_size);

        // Step in each space direction
        for (int j = 0; j < space_dim; j++)
        {
            // Random step in j direction between -delta and + delta
            double dj = ((drand48() * 2) - 1) * delta;

            // NOTE: In periodic boundaries conditions there is no need to check for boundaries
            pos_array[c(i, j)] += dj;
            steps_save[j] = dj;
        }

        // Compute the new energy
        double new_energy = compute_i_lennar_jones_potential(i, pos_array, charge_array, n_particles, space_dim, box_size);
        double dE = new_energy - old_energy;

        // Metropolis acceptance criterion
        int accepted = 0;
        double alpha = fmin(1, exp(-dE / temperature));
        accepted = drand48() <= alpha;

        if (accepted)
        {
            (*accepted_counter)++;

            // Update energy
            energy += dE;

            // Bring particles back to first cell
            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] = pb_wrap_position(pos_array[c(i, j)], box_size);
            }
        }
        else
        {
            // Restore positions
            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] -= steps_save[j];
            }
        }
    }

    return energy;
}

int main(int argc, char const *argv[])
{
    //---------------------------
    // SIMULATION PARAMETERS
    //---------------------------
    const int lattice_type = 4;
    const int n_cell_per_row = 5;
    const double density = 0.85;
    //---------------------------
    const int space_dimension = 3;
    const int seed = 42;
    const int n_metropolis_step = 100000;
    const double space_step = 0.1;
    //---------------------------
    // DERIVATE PARAMETERS
    //---------------------------
    const int n_particles = pow(n_cell_per_row, 3) * lattice_type;
    const int total_vel_pos_array_size = n_particles * space_dimension;
    const double box_size = pow(n_particles / density, 1.0 / space_dimension);
    const int PRINT_INTERVAL = n_metropolis_step / 100; // Update progress bar every 1%
    //---------------------------

    printf("N_particles = %d\n", n_particles);
    printf("Box_size = %lf\n", box_size);

    srand48(seed);

    // g(r) CONSTANTS --------------------
    const double max_radius = box_size; // After box_size/2 one can see artifact emerges from the conflict between the spherical nature of g(r) and the cubic cell
    const double N_bins = 600;
    const double bin_interval = max_radius / N_bins;
    int n_radial_distributions_performed = 0; // Keep track of how many times the g(r) is computed and summed to the bins, used to normalize

    // Array of bins used to compute g(r)
    double *bin_counting_array = (double *)calloc(sizeof(double), N_bins); // Array di bin per tenere traccia della densità radiale
    if (!bin_counting_array)
    {
        perror("Error allocating counting array");
        exit(EXIT_FAILURE);
    }

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

    double *charge_array = (double *)malloc(n_particles * sizeof(double));
    if (charge_array == NULL)
        exit(EXIT_FAILURE);

    double *mass_array = (double *)malloc(n_particles * sizeof(double));
    if (mass_array == NULL)
        exit(EXIT_FAILURE);

    FILE *start_position_file = fopen("./output/start_position_file.csv", "w");
    if (start_position_file == NULL)
    {
        printf("Check for output folder");
        exit(EXIT_FAILURE);
    }

    FILE *end_position_file = fopen("./output/end_position_file.csv", "w");
    if (end_position_file == NULL)
        exit(EXIT_FAILURE);

    FILE *radial_distribution_file = fopen("./output/radial_distribution.csv", "w");
    if (radial_distribution_file == NULL)
        exit(EXIT_FAILURE);

    FILE *energy_file = fopen("./output/energy.csv", "w");
    if (energy_file == NULL)
        exit(EXIT_FAILURE);

    setvbuf(energy_file, NULL, _IOFBF, 1 << 20); // 1 MB buffer

    //-------------------------------------------

    init_system_lattice(pos_array, charge_array, mass_array, n_particles, box_size, lattice_type, n_cell_per_row);
    // init_system(pos_array, charge_array, mass_array, n_particles, space_dimension, box_size);

    // Compute energy in the starting position
    double energy = pb_compute_total_energy(pos_array, charge_array, n_particles, space_dimension, box_size);

    save_particle_state(start_position_file, pos_array, charge_array, n_particles, space_dimension);
    printf("Start energy: %f\n", energy);

    clock_t begin_time = clock(); // Save starting time, used for ETA in progress bar

    long accepted_steps = 0; // Keep track of how many metropolis steps get accepted

    // INCREASE T METHOD -------------------------------
    const double temp_min = 0.45;
    const double temp_max = 0.90;
    const int N_temperatures = 10;
    const double dT = (temp_max - temp_min) / N_temperatures;

    const int N_step_thermalization = 1000;
    const int N_step_data = 10000;

    double *energy_array = (double *)malloc(sizeof(double) * N_step_data);
    if (energy_array == NULL)
        exit(EXIT_FAILURE);

    // Array to keep track of all specific_heat varing temperature
    double *specific_heat_array = (double *)malloc(sizeof(double) * N_temperatures);
    if (specific_heat_array == NULL)
        exit(EXIT_FAILURE);

    FILE *specific_heat_file = fopen("./output/specific_heat.csv", "w");
    if (specific_heat_file == NULL)
        exit(EXIT_FAILURE);

    for (size_t i = 0; i < N_temperatures; i++)
    {
        double temperature = temp_min + i * dT;
        printf("T = %lf\n", temperature);

        // Termalizations steps plus data steps
        for (size_t j = 0; j < N_step_thermalization + N_step_data; j++)
        {
            // Progress bar
            if (j % ((N_step_thermalization + N_step_data) / 100) == 0)
            {
                print_progress(j, N_step_thermalization + N_step_data, begin_time);
                printf(" %d/%d", (int)i + 1, N_temperatures);
                fflush(stdout);
            }

            // Perform metropolis step
            energy = pb_metropolis_step_one_particle(energy, pos_array, charge_array, space_step, temperature, n_particles, space_dimension, &accepted_steps, box_size);

            fprintf(energy_file, "%lf\n", energy);

            // Start keep track of energy for statistics only after thermalization
            if (j >= N_step_thermalization)
            {
                energy_array[j - N_step_thermalization] = energy;
            }
        }

        fprintf(specific_heat_file, "%lf;%lf\n", temperature, array_var(energy_array, N_step_data) / (temperature * temperature) + 3. / 2.0 * n_particles);

        // Flush data in order to avoid data loss in the middle of the simulation
        fflush(specific_heat_file);
        fflush(energy_file);

        printf("\r\033[2K"); // Erase cmd current line
    }

    // END INCREASE T METHOD -------------------------------

    //
    // for (size_t i = 0; i < n_metropolis_step; i++)
    //{
    //    // energy = pb_metropolis_step_full_system(energy, pos_array, charge_array, space_step, temperature, n_particles, space_dimension, &accepted_steps, box_size);
    //    energy = pb_metropolis_step_one_particle(energy, pos_array, charge_array, space_step, temperature, n_particles, space_dimension, &accepted_steps, box_size);
    //
    //    if (i % PRINT_INTERVAL == 0)
    //    {
    //        // Progress Bar
    //        print_progress(i, n_metropolis_step, begin);
    //        fflush(energy_file);
    //
    //        // Trashold for termalization
    //        if (i > 1000)
    //        {
    //            n_radial_distributions_performed++;
    //            radial_distribution(pos_array, n_particles, space_dimension, box_size, bin_counting_array, N_bins, bin_interval);
    //        }
    //    }
    //
    //    if (i % 10 == 0)
    //    {
    //        fprintf(energy_file, "%lf\n", energy);
    //    }
    //}
    //
    // Clear terminal
    printf("\r\033[2K");

    //----------------------------------
    clock_t end_time = clock();
    double time_spent = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
    printf("Total time: %.0lf ms\n", time_spent * 1000);

    printf("Accepted step: %ld\n", accepted_steps);
    printf("End energy: %f\n", energy);

    save_particle_state(end_position_file, pos_array, charge_array, n_particles, space_dimension);

    // Normalization of bins in g(r)
    for (size_t i = 0; i < N_bins; i++)
    {
        bin_counting_array[i] = (double)bin_counting_array[i] / (density * n_particles * n_radial_distributions_performed * (4 * PI / 3) * (pow((i + 1) * bin_interval, 3) - pow(i * bin_interval, 3)));
    }

    for (size_t j = 0; j < N_bins; j++)
    {
        fprintf(radial_distribution_file, "%f;%f\n", j * bin_interval / (box_size / 2), bin_counting_array[j]);
    }

    free(pos_array);
    free(vel_array);
    free(mass_array);
    free(charge_array);
    free(specific_heat_array);
    fclose(start_position_file);
    fclose(end_position_file);
    fclose(energy_file);
    fclose(specific_heat_file);
    fclose(radial_distribution_file);

    return 0;
}
