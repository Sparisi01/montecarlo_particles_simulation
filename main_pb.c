/**
 * @file main_pb.c
 * @brief
 *
 * @details
 *
 * @author Parisi Simone and Lenz Patrick
 * @date 17 December 2025
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>

#include "src/constants.c"
#include "src/checkpoints_handler.c"
#include "src/arrays_stat_operations.c"
#include "src/lennar_jones.c"
#include "src/periodic_boundaries.c"
#include "src/ewald.c"
#include "src/progress_bar.c"
#include "src/verlet_list.c"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIX(a, b) ((a) < (b) ? (a) : (b))

enum SIMULATION_TYPE
{
    SINGLE_T,
    INCREASING_T,
    NONE
};

const int COULOMB_INTERACTION_ON = 0;

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

/**
 * @brief Save particles position and charge state in a csv file, easy to read in python for data analysis.
 */
void save_particle_state_csv(const char *filename,
                             double *pos_array,
                             double *charge_array,
                             int n_particles,
                             int space_dim)
{
    FILE *f = fopen(filename, "w");
    if (!f)
    {
        perror("Csv save failed");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dim; j++)
        {
            fprintf(f, "%f;", pos_array[c(i, j)]);
        }
        fprintf(f, "%f\n", charge_array[i]);
    }

    fclose(f);
}

/**
 * @brief Given the system of N particle in D dimensional space compute the total interaction energy under periodic boudary conditions.
 * By default only Lennar Jones potential is used. If at least 2 particles have a charge different from 0 the coulomb potential is turned on.
 *
 * @note Coulomb potential is computed used Ewald Summation and it requires a space dimension of 3.
 */
double pb_compute_total_energy(const double *pos_array,
                               const double *charge_array,
                               int n_particles,
                               int space_dim,
                               double box_size,
                               double epsilon,
                               double sigma)
{
    double total_energy = 0;
    total_energy += pb_total_lennar_jones_energy(pos_array, charge_array, n_particles, space_dim, box_size, epsilon, sigma);

    if (COULOMB_INTERACTION_ON)
    {
        if (space_dim != 3)
        {
            fprintf(stderr, "Error: Ewald Summation require a space dimension of 3\n");
        }
        total_energy += ewd_total_energy(pos_array, charge_array, n_particles, box_size);
    }

    return total_energy;
}

/**
 * @brief Same as "pb_compute_total_energy" but using VerletList to decrease computations efford. See "pb_compute_total_energy"
 * for more on comments.
 */
double pb_verlet_compute_total_energy(const double *pos_array,
                                      const double *charge_array,
                                      const IndexesList_t *verlet_list,
                                      int n_particles,
                                      int space_dim,
                                      double box_size,
                                      double epsilon,
                                      double sigma)
{
    double total_energy = 0;
    total_energy += pb_verlet_tot_lennar_jones_energy(pos_array, charge_array, verlet_list, n_particles, space_dim, box_size, epsilon, sigma);

    if (COULOMB_INTERACTION_ON)
    {
        if (space_dim != 3)
        {
            fprintf(stderr, "Error: Ewald Summation require a space dimension of 3\n");
        }
        total_energy += ewd_verlet_total_energy(pos_array, charge_array, verlet_list, n_particles, box_size);
    }

    return total_energy;
}

// Initialize positions in a cubic lattice
void init_system_lattice(double *pos_array,
                         double *charge_array,
                         double *mass_array,
                         int n_particles,
                         double box_size,
                         int lattice_type,
                         int n_cell_per_row)
{
    struct Vec3
    {
        double x;
        double y;
        double z;
    } typedef Vec3;

    // lattice_type 1 CC, 2 BCC, 4 FCC

    Vec3 lattice_positions[lattice_type];
    int basis_parity[lattice_type];

    switch (lattice_type)
    {
    case 1: // CC
        lattice_positions[0] = (Vec3){0, 0, 0};
        basis_parity[0] = 0;
        basis_parity[1] = 1;
        basis_parity[2] = 1;
        basis_parity[3] = 1;
        break;
    case 2: // BCC
        lattice_positions[0] = (Vec3){0, 0, 0};
        lattice_positions[1] = (Vec3){0.5, 0.5, 0.5};
        basis_parity[0] = 0;
        basis_parity[1] = 1;
        break;
    case 4: // FCC
        lattice_positions[0] = (Vec3){0, 0, 0};
        lattice_positions[1] = (Vec3){0, 0.5, 0.5};
        lattice_positions[2] = (Vec3){0.5, 0.5, 0};
        lattice_positions[3] = (Vec3){0.5, 0, 0.5};
        basis_parity[0] = 0;
        basis_parity[1] = 1;
        basis_parity[2] = 1;
        basis_parity[3] = 1;
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

                    int parity = (n_cell_x + n_cell_y + n_cell_z + basis_parity[w]) % 2;
                    charge_array[n_particle_placed] = (parity == 0) ? 1.0 : -1.0;

                    mass_array[n_particle_placed] = 1;
                    n_particle_placed++;
                }
            }
        }
    }
}

void init_system(double *pos_array,
                 double *charge_array,
                 double *mass_array,
                 int n_particles,
                 int space_dimension,
                 double box_size)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = 0; j < space_dimension; j++)
        {
            pos_array[c(i, j)] = drand48() * box_size; // Uniform position distribution inside the square box
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

// This function perform an update of all the system at once, hard to get good acceptance rate for stable configurations.
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
    double new_energy = pb_compute_total_energy(pos_array, charge_array, n_particles, space_dim, box_size, 1, 1);
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

// This function perform an update of all the system on one particle at the time, easier to get good acceptance rate but
// it has highter temporal correlation
double verlet_pb_metropolis_step_one_particle(double energy,
                                              double *pos_array,
                                              const double *charge_array,
                                              const IndexesList_t *vl,
                                              double delta,
                                              double temperature,
                                              int n_particles,
                                              int space_dim,
                                              long *accepted_counter,
                                              double box_size,
                                              double epsilon,
                                              double sigma)
{

    /** In order to remove the bias caused by updating all the particle using always the
     * same order we define a permutation array of indexes. At each metropolis step
     * a new random permutation is used, ensuring that all the particles are updated but in
     * different order. Based on literature this ensured a drop in time correlation.
     */
    static int *perm = NULL;
    if (perm == NULL)
    {
        perm = (int *)malloc(n_particles * sizeof(int));
        for (int i = 0; i < n_particles; i++)
        {
            perm[i] = i;
        }
    }

    // NOTE: I have tested with and without this line of code and we move from an
    // autocorrelation af 12 steps with to 145 without

    for (int i = n_particles - 1; i > 0; i--)
    {
        // Swap two random indexes in the sub array 0,1,...,i
        int j = (int)(drand48() * (i + 1));
        int tmp = perm[i];
        perm[i] = perm[j];
        perm[j] = tmp;
    }

    static double *old_position = NULL;
    if (old_position == NULL)
    {
        old_position = (double *)malloc(sizeof(double) * n_particles * space_dim);
    }

    memcpy(old_position, pos_array, n_particles * space_dim);

    // Keep track of space shifts for the current particle
    double steps_save[space_dim];

    for (size_t p = 0; p < n_particles; p++)
    {

        int i = perm[p];

        double old_energy = pb_verlet_i_lennar_jones_potential(i, pos_array, charge_array, vl, n_particles, space_dim, box_size, epsilon, sigma);
        if (COULOMB_INTERACTION_ON)
            old_energy += ewd_verlet_i_short_energy(i, pos_array, charge_array, vl, n_particles, box_size);

        // Random step in j direction between -delta and + delta
        for (int j = 0; j < space_dim; j++)
        {
            double dj = ((drand48() * 2) - 1) * delta;

            pos_array[c(i, j)] += dj;
            steps_save[j] = dj;
        }

        double new_energy = pb_verlet_i_lennar_jones_potential(i, pos_array, charge_array, vl, n_particles, space_dim, box_size, epsilon, sigma);
        if (COULOMB_INTERACTION_ON)
            new_energy += ewd_verlet_i_short_energy(i, pos_array, charge_array, vl, n_particles, box_size);

        double dE = new_energy - old_energy;

        if (COULOMB_INTERACTION_ON)
            dE += ewd_delta_long_energy(i, pos_array, old_position, charge_array, box_size);

        /** Metropolis acceptance criterion
         *
         * @brief If the step is accepted update the current total energy using dE
         * and apply periodic boundary condition bringing all the particle back
         * to the first cubic cell.
         *
         * If the step is refused, undo the space step using the steps dj
         * saved in steps_save.
         */
        int accepted = 0;
        double alpha = fmin(1, exp(-dE / temperature));
        accepted = drand48() <= alpha;

        if (accepted)
        {
            (*accepted_counter)++;

            energy += dE;

            if (COULOMB_INTERACTION_ON)
                ewd_update_S_K(i, pos_array, old_position, charge_array, box_size);

            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] = pb_wrap_position(pos_array[c(i, j)], box_size);
                old_position[c(i, j)] = pos_array[c(i, j)];
            }
        }
        else
        {
            for (int j = 0; j < space_dim; j++)
            {
                pos_array[c(i, j)] -= steps_save[j];
            }
        }
    }

    return energy;
}

void print_simulation_information(const int n_particles,
                                  const double box_size,
                                  IndexesList_t *verlet_list,
                                  double energy,
                                  double *pos_array,
                                  double *charge_array,
                                  const int space_dimension,
                                  double density)
{
    printf("\n");
    printf("=====================================\n");
    printf("        STARTING SIMULATION\n");
    printf("=====================================\n");

    printf(" N_particles : %d\n", n_particles);
    printf(" Box_size    : %.3f\n", box_size);
    printf(" Density     : %.3f\n", density);

    printf("-------------------------------------\n");
    printf(" Max verlet count: %d\n", get_max_verlet_count(verlet_list, n_particles));
    printf("-------------------------------------\n");

    printf(" Start energy      : %.6e\n", energy);

    double no_verlet_energy = pb_compute_total_energy(pos_array, charge_array, n_particles, space_dimension, box_size, 1, 1);
    assert(fabs(energy - no_verlet_energy) / fabs(energy) < 1e-3);

    printf("-------------------------------------\n");
}

// ------------------------------------------------------------------

int main(int argc, char const *argv[])
{
    if (COULOMB_INTERACTION_ON)
        printf("COULOMB INTERACTION ON\n");
    else
        printf("COULOMB INTERACTION OFF\n");

    /**
     * SIMULATION PARAMETERS
     * Those are the 3 parameters that actually control the simulation.
     * lattice_type & n_cell_per_row define the number of particles,
     * density and number of particles define the box size.
     */
    const int lattice_type = 4;   // Lattice type 1 CC, 2 BCC, 4 FCC
    const int n_cell_per_row = 6; // Number of lattice cell per row
    const double density = 1;

    const int space_dimension = 3; // 1D - 2D - 3D - ... - nD

    const int seed = 42;
    srand48(seed);

    const int n_particles = pow(n_cell_per_row, 3) * lattice_type;
    const double box_size = pow(n_particles / density, 1.0 / space_dimension);

    const double space_step = 0.1; // Max space step in metropolis update, should be << box_size

    if (space_step > box_size * 0.05)
    {
        fprintf(stderr, "WARNING: space_step could be too hight for the current box_size.\n");
        fprintf(stderr, "Space_step = %lf, Box_size = %lf\n", space_step, box_size);
    }

    const double lennar_jones_epsilon = 1;
    const double lennar_jones_sigma = 1;

    /**
     * Positions and velocities are store in a flattered array of size (N * SPACE_DIM), where
     * the j component of the i particle is stored at index (SPACE_DIM * i + j).
     * In order to retrieve the correct index for component j of particle i use the macro "c(i,j)" defined
     * in constants.c.
     */
    const int total_vel_pos_array_size = n_particles * space_dimension;

    double *pos_array = (double *)malloc(total_vel_pos_array_size * sizeof(double));
    if (pos_array == NULL)
        exit(EXIT_FAILURE);

    double *charge_array = (double *)malloc(n_particles * sizeof(double));
    if (charge_array == NULL)
        exit(EXIT_FAILURE);

    double *mass_array = (double *)malloc(n_particles * sizeof(double));
    if (mass_array == NULL)
        exit(EXIT_FAILURE);

    FILE *radial_distribution_file = fopen("./output/radial_distribution.csv", "w");
    if (radial_distribution_file == NULL)
        exit(EXIT_FAILURE);

    FILE *energy_file = fopen("./output/energy.csv", "w");
    if (energy_file == NULL)
        exit(EXIT_FAILURE);

    setvbuf(energy_file, NULL, _IOFBF, 1 << 20); // 1 MB buffer, used to increase performance

    /**
     * And estimate of g(r) is computed using binning,
     * the following parameters define the binning number and size
     */
    const double max_radius = box_size / 2; // After box_size/2 one can see artifact emerges from the conflict between the spherical nature of g(r) and the cubic cell
    const int N_bins = 600;
    const double bin_interval = max_radius / N_bins;
    int n_radial_distributions_performed = 0; // Keep track of how many times the g(r) is computed and summed to the bins, used to normalize

    double *bin_counting_array = (double *)calloc(N_bins, sizeof(double)); // Bins array
    if (!bin_counting_array)
    {
        perror("Error allocating counting array");
        exit(EXIT_FAILURE);
    }

    long metropolis_accepted_steps = 0;
    double energy = 0;

    const int restart_from_checkpoint = 1;
    if (restart_from_checkpoint)
    {
        const char *check_point_file_name = "./state_saves_binaries/checkpoint_T1.500.bin";
        load_checkpoint_binary(check_point_file_name, pos_array, n_particles, space_dimension, &energy);
    }
    else
    {
        init_system_lattice(pos_array, charge_array, mass_array, n_particles, box_size, lattice_type, n_cell_per_row);
        // init_system(pos_array, charge_array, mass_array, n_particles, space_dimension, box_size);
    }

    save_particle_state_csv("./output/start_position_file.csv", pos_array, charge_array, n_particles, space_dimension);

    /** |------VERLET LIST------|
     *
     * @brief The energy computation is O(N^2). O(N) for each particle and the total is N * O(N).
     * This is not feasable for large N, so we use a techique called verlet List, which we will briefly explain here.
     *
     * The short range potentials like Lennar Jones and real parte of Ewald Summation can be truncated without loss
     * of accuracy, each interaction above r = r_c are not computed. Using this fact the idea is the following:
     * for each particle "i" we define an array of indexes, in this array are contained the indexes of all the particle
     * that are no more than r_c + r_skin away from "i". When the energy has to be computed we can loop over these lists
     * instead of all the particles.
     * r_skin is used to ensure that we can wait a certain number of metropolis step before rebuild the lists (Remember that
     * building the list is actually O(N^2), if we have to do it at every step the would be no increase in performance).
     *
     * The alghorithm using verlet is the following:
     *
     * 1) given the particle array build N lists, one for each particle.
     * 2) compute the energy looping over the lists.
     * 3) perform a metropolis step and then check if any of the particle has move more than r_skin (this checking is O(N)),
     * if it's the case we need to update the list, if not we can perform poit 3) again.
     *
     * The hard part is to balance the frequency of list rebuild and number of particle per list. One in theory wants low
     * rebuild frequency and low number of particle per list. This is the case for low density states where the particle are
     * far apart they take a lot of metropolis steps to actually move near another particle.
     *
     * NOTE: is there a way to compute the optimal skin based on density, N and step_size?
     * I think one can do that by approximating the metropolis step as a random walk and using the fact that
     * we how which std it has, sqrt(N).
     */
    const double verlet_max_neightbor_distance = 2.5 * SIGMA; // Same used in Lennar Jones
    const double skin = 0.5 * verlet_max_neightbor_distance;

    double *old_pos_array = (double *)malloc(total_vel_pos_array_size * sizeof(double));
    if (old_pos_array == NULL)
    {
        exit(EXIT_FAILURE);
    }

    IndexesList_t *verlet_list = (IndexesList_t *)malloc(sizeof(IndexesList_t) * n_particles);
    if (verlet_list == NULL)
    {
        exit(EXIT_FAILURE);
    }

    // Build verlet list
    verlet_pb_build_list(pos_array, old_pos_array, verlet_list, n_particles, space_dimension, box_size, verlet_max_neightbor_distance, skin);

    // Optimize Ewald Summation Parameters
    if (COULOMB_INTERACTION_ON)
    {
        optimizeParameter(1, box_size, charge_array, n_particles);
        if (REAL_CUTOFF > verlet_max_neightbor_distance)
            LOG_WARNING("REAL_CUTOFF greater than verlet_max_neightbor_distance.\n"
                        "There is the possibility that some ral space coulomb interactions may be excluded by the verlet list.\n"
                        "Consider increasing verlet_max_neightbor_distance or decrease the REAL_CUTOFF by optimizing ewald with using a larger error.");
    }

    // Initialization of S(K) vector in the starting position
    ewd_init_S_K(pos_array, charge_array, n_particles, box_size);

    if (!restart_from_checkpoint)
    {
        energy = pb_verlet_compute_total_energy(pos_array, charge_array, verlet_list, n_particles, space_dimension, box_size, lennar_jones_epsilon, lennar_jones_sigma);
    }

    // Print a lot of informations about the simulation in order to spot possible errors at the start of simulation
    print_simulation_information(n_particles, box_size, verlet_list, energy, pos_array, charge_array, space_dimension, density);

    // Choose type of simulation
    enum SIMULATION_TYPE simulation_type = INCREASING_T;

    switch (simulation_type)
    {
    case INCREASING_T:
        goto INCREASE_TEMPERATURE_SIMULATION;
        break;

    case SINGLE_T:
        goto SINGLE_TEMPERATURE_SIMULATION;
        break;

    case NONE:
        goto FREE_SECTION;
        break;
    }

INCREASE_TEMPERATURE_SIMULATION:

    /** |---- INCREASE TEMPERATURE SIMULATION -----|
     * Start from a system configuration and define a starting temperature T, a temperature step dT
     * and a number of temperature step N (T_max = Tmin + NdT)
     *
     * 1) Set the min temperature.
     *
     * 2) The system is evolved throught metropolis for N_step_thermalization until thermalization.
     *
     * 3) After the termalization the system is evolved for N_step_data steps. During this time
     * the thermodinamics varables are stored (Energy). At the end compute statistics on the current
     * energy array (Cv).
     *
     * 4) Save the current system state ina  binary file
     *
     * 5) Increase T by dT and repeat
     *
     *
     * NOTE: The total number of steps is (N_step_thermalization + N_step_data) * N_temperatures
     */

    const double temperature_min = 1.50;
    const double temperature_max = 1.60;
    const int N_temperatures = 5;
    const double dT = (temperature_max - temperature_min) / N_temperatures;

    const int N_step_thermalization = 5000;
    const int N_step_data = 50000;
    const int N_step_tot = N_step_thermalization + N_step_data;

    double *energy_array = (double *)malloc(sizeof(double) * N_step_data);
    if (energy_array == NULL)
        exit(EXIT_FAILURE);

    double *specific_heat_array = (double *)malloc(sizeof(double) * N_temperatures);
    if (specific_heat_array == NULL)
        exit(EXIT_FAILURE);

    FILE *specific_heat_file = fopen("./output/specific_heat.csv", "w");
    if (specific_heat_file == NULL)
        exit(EXIT_FAILURE);

    for (size_t i = 0; i < N_temperatures; i++)
    {
        clock_t begin_time = clock(); // Save starting time, used for ETA in progress bar

        double temperature = temperature_min + i * dT;
        printf("T = %lf\n", temperature);

        // Termalizations steps plus data steps
        for (size_t j = 0; j < N_step_tot; j++)
        {
            // Progress bar
            if (j % (N_step_tot / 100) == 0)
            {
                print_progress(j, N_step_tot, begin_time);
                printf(" %d/%d", (int)i + 1, N_temperatures);
                fflush(stdout);
            }

            energy = verlet_pb_metropolis_step_one_particle(energy, pos_array, charge_array, verlet_list, space_step, temperature, n_particles, space_dimension, &metropolis_accepted_steps, box_size, lennar_jones_epsilon, lennar_jones_sigma);

            // TODO: add verlet list update and usage

            fprintf(energy_file, "%lf\n", energy);

            // Start keep track of energy for statistics only after thermalization
            if (j >= N_step_thermalization)
            {
                energy_array[j - N_step_thermalization] = energy;
            }
        }

        fprintf(specific_heat_file, "%lf;%lf\n", temperature, array_var(energy_array, N_step_data) / (temperature * temperature));

        // Flush data in order to avoid data loss in the middle of the simulation
        fflush(specific_heat_file);
        fflush(energy_file);

        printf("\r\033[2K"); // Erase cmd current line

        // Save check point at current T
        char filename[128];
        sprintf(filename, "./state_saves_binaries/checkpoint_T%.3f.bin", temperature);
        save_checkpoint_binary(filename, pos_array, n_particles, space_dimension, energy);
    }

    free(energy_array);
    free(specific_heat_array);
    fclose(specific_heat_file);

    goto FREE_SECTION;

    // |---- END INCREASE TEMPERATURE SIMULATION -----|

SINGLE_TEMPERATURE_SIMULATION:

    /** |------ SINGLE TEMPERATURE SIMULATION ------|
     *
     * Starting from a start position for all the particle N metropolis steps are performed
     * keeping the temperature fixed. Usefull to study the equilibrium radial distribution and
     * equilibrium thermodinamic quantities.
     *
     * Note that all the statistics has to be performed AFTER reaching equilibrum and keeping in mind that
     * there is time correlation between metropolis steps. This is not a huge problem for g(r) but it is
     * when we compute the std of thermodinamic quantities.
     *
     * Following test has been succesfully performed (all starting from FCC lattice with only lennar jones potential):
     *
     * - (density = 0.1, T = 1.1) gas behaviour can be observed. g(r) has a strong peak around 1 (the minimum of lennar jones potential) and
     * for hight distances the density is constant.
     * - (density = 0.7, T = 1.1) liquid behaviour can be observed. g(r) has a few shallow peaks.
     * - (density = 1.3, T = 1.1) solid behaviour can be observed. g(r) has a lot of strong peaks.
     *
     * Only LJ varing T fixes rho
     * - (density = 1, N = 864, T = 1)  Solid
     * - (density = 1, N = 864, T = 1.5)  Solid
     * - (density = 1, N = 864, T = 2)  Liquid
     */

    const int N_data_steps = 10000;
    const int N_thermalization_steps = 5000;
    const int N_metropolis_steps = N_thermalization_steps + N_data_steps;
    double temperature = 2;

    printf(" Temperature : %.3f\n", temperature);
    printf("-------------------------------------\n");

    // Keep track of how frequntly the verlet lists are updated
    int counting_verlet = 0;
    int min_counting_verlet = INT32_MAX;

    clock_t begin_time = clock();
    for (size_t i = 0; i < N_metropolis_steps; i++)
    {
        counting_verlet++;

        // energy = pb_metropolis_step_full_system(energy, pos_array, charge_array, space_step, temperature, n_particles, space_dimension, &accepted_steps, box_size);
        energy = verlet_pb_metropolis_step_one_particle(energy, pos_array, charge_array, verlet_list, space_step, temperature, n_particles, space_dimension, &metropolis_accepted_steps, box_size, lennar_jones_epsilon, lennar_jones_sigma);

        // Check if the verlet list need to be rebuild
        if (verlet_pb_needs_rebuild(pos_array, old_pos_array, n_particles, space_dimension, box_size, skin))
        {
            verlet_pb_build_list(pos_array, old_pos_array, verlet_list, n_particles, space_dimension, box_size, verlet_max_neightbor_distance, skin);
            if (counting_verlet < min_counting_verlet)
            {
                min_counting_verlet = counting_verlet;
            }
            counting_verlet = 0;
        }

        if (i % (N_metropolis_steps / 100) == 0)
        {
            // Progress Bar
            print_progress(i, N_metropolis_steps, begin_time);
            printf(" - Min number of step before update: %d", min_counting_verlet);

            fflush(stdout);
            fflush(energy_file);
        }

        // Trashold for termalization
        if ((i > N_thermalization_steps) & (i % 100 == 0))
        {
            n_radial_distributions_performed++;
            radial_distribution(pos_array, n_particles, space_dimension, box_size, bin_counting_array, N_bins, bin_interval);
        }

        if (i > 0)
        {
            fprintf(energy_file, "%lf\n", energy);
        }
    }

    // Clear terminal and print end simulation info
    printf("\r\033[2K");
    printf("Min number of step before update: %d\n", min_counting_verlet);
    printf("Accepted step: %ld\n", metropolis_accepted_steps);
    printf("End energy: %f\n", energy);
    printf("=====================================\n");

    save_particle_state_csv("./output/end_position_file.csv", pos_array, charge_array, n_particles, space_dimension);

    // Normalization of bins in g(r)
    for (size_t i = 0; i < N_bins; i++)
    {
        bin_counting_array[i] = (double)bin_counting_array[i] / (density * n_particles * n_radial_distributions_performed * (4 * PI / 3) * (pow((i + 1) * bin_interval, 3) - pow(i * bin_interval, 3)));
    }

    for (size_t j = 0; j < N_bins; j++)
    {
        double r = (j + 0.5) * bin_interval;
        fprintf(radial_distribution_file, "%f;%f\n", r, bin_counting_array[j]);
    }

    // |------ END SINGLE TEMPERATURE SIMULATION ------|

FREE_SECTION:

    free(pos_array);
    free(mass_array);
    free(charge_array);
    free(verlet_list);
    free(old_pos_array);

    fclose(energy_file);
    fclose(radial_distribution_file);

    return 0;
}
