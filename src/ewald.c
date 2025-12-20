/**
 * @file ewald.c
 * @brief framework to compute coulomb like energy in a 3D periodic boundary condition eviroment.
 * Both with and without verlet list implementations are given.
 *
 * The following files are required:
 * "constants.c"
 * "periodic_boundaries.c"
 * "verlet_list.c"
 *
 * @cite TODO:
 *
 * @author Parisi Simone
 * @date 17 December 2025
 */

#ifndef EWALD_H
#define EWALD_H

#include <math.h>
#include <inttypes.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <string.h>

#include "constants.c"
#include "periodic_boundaries.c"
#include "verlet_list.c"
#include "logger.c"

#define EWD_EPSILON 1e-8

static double RECIPROCAL_RANGE;
static double REAL_CUTOFF;
static double ALPHA;
static int OPTIMIZED = 0;

static double complex *S_K = NULL;

const char CUTOFF_WARNING_MESSAGE[] =
    "REAL_CUTOFF exceeds box_size/2, violating the first-image convention "
    "assumed by \"ewd_verlet_i_short_energy\". "
    "Computed energies may be incorrect. "
    "Consider increasing the \"error\" tolerance or manually tuning "
    "ALPHA, REAL_CUTOFF, and RECIPROCAL_RANGE.\n";

/**
 * @brief The following three functions perform an optimization of Ewald parameters given
 *
 * - ABSOLUTE error on the energy
 * - box_size
 * - charge_array
 * - n_particles
 *
 * these optimized parameters, used in combination with verlet list, bring the complexity
 * of one ewd_total_energy from O(N) to O(N^3/2).
 *
 * @cite TODO: cite the paper for this part
 */

static double errorsDifference(double error, double s, double Q2, double cell_length)
{
    return exp(-(s * s)) / (pow(s, 3. / 2)) * Q2 * sqrt((2 + PI) / (2 * PI * ALPHA * pow(cell_length, 3))) - error;
};

static double findSbybisection(double a, double b, double error, double Q2, double cell_length, double precision)
{
    double max = errorsDifference(error, a, Q2, cell_length);
    double min = errorsDifference(error, b, Q2, cell_length);

    if (!(max > 0 && min < 0))
    {
        LOG_ERROR("Max e min non rispettano parametri bisezione");
    };
    double c = 0;
    int root_find = 0;
    while (!root_find)
    {
        c = (a + b) / 2;
        double mid_value = errorsDifference(error, c, Q2, cell_length);
        if (fabs(mid_value) < precision)
        {
            root_find = 1;
        }
        else
        {
            if (mid_value < 0)
                b = c;
            if (mid_value > 0)
                a = c;
        }
    }
    return c;
}

void optimizeParameter(double error, double box_size, const double *charge_array, int n_particles)
{

    double Q2 = 0;
    for (size_t i = 0; i < n_particles; i++)
    {
        Q2 += charge_array[i] * charge_array[i];
    }

    if (Q2 == 0)
    {
        fprintf(stderr, "\033[1;31mERROR: All charges are set to 0. Optimization cannot proceed with zero charges.\n");
        fprintf(stderr, "Possible solutions:\n");
        fprintf(stderr, "- Disable Coulomb interaction if it's not needed.\n");
        fprintf(stderr, "- Modify the 'charge_array' to use non-zero charge values.\n");
        fprintf(stderr, "- Ensure the charge values are correctly initialized.\n");
        fprintf(stderr, "Please check your input parameters and try again.\n");
        fprintf(stderr, "\033[0m"); // Resets color to default
    }

    // These 3 parameters actually depends on the machine you are running these, but the result for ALPHA is not hightly sensible on TAU_RAPP
    const double TAU_S = 3.6; // Time needed to perform one short range energy evaluation (per particle)
    const double TAU_L = 1.0; // Time needed to perform one long range energy evaluation (per particle)
    const double TAU_RAPP = (TAU_S / TAU_L);

    ALPHA = pow((TAU_RAPP * pow(PI, 3) / pow(box_size, 6) * n_particles), 1. / 6);

    double s = findSbybisection(1e-3, 1e3, error, Q2, box_size, 1e-10);
    REAL_CUTOFF = s / ALPHA;

    double kc = 2 * s * ALPHA;
    RECIPROCAL_RANGE = ceil(kc / (2 * PI / box_size));

    printf("OPTIMIZED PARAMETERS: ALPHA = %.5E, R_C = %.5E, N_C_K = %.5E\n", ALPHA, REAL_CUTOFF, RECIPROCAL_RANGE);

    if (REAL_CUTOFF > box_size / 2)
    {
        LOG_WARNING("%s", CUTOFF_WARNING_MESSAGE);
    }

    OPTIMIZED = 1;
}

// Much simple parameter oprimization
static inline void ewd_alpha_by_precision(double precision)
{
    ALPHA = sqrt(-log(precision));
}

/**
 * @brief Compute the self-energy of a system of charged particles.
 *
 * Calculates the self-energy contribution of the Coulomb interaction for each particle
 * in the system. This term represents the energy of a particle due to its own charge
 * distribution in the system. The self-energy is computed by summing the squares of
 * the charges of all particles, applying the screening parameter `alpha`, and normalizing
 * by a factor of `1 / √π`.
 */
double ewd_self_energy(const double *charge_array, int n_particles)
{
    double sum = 0;
    for (size_t i = 0; i < n_particles; i++)
    {
        sum += charge_array[i] * charge_array[i];
    }
    return sum * ALPHA / SQRT_PI;
}

/**
 * @brief Compute the real-space Coulomb energy of a particle i in a system of charged particles.
 *
 * Calculates the short-range Coulomb interactions between a particle i and all others
 * in the system.
 * The interactions are truncated at a distance defined by the cutoff (_CUTOFF) to limit
 * long-range interactions, and the Coulomb potential is screened by the complementary
 * error function (`erfc`) with a decay parameter `alpha`.
 *
 * @note Complexity O(N)
 */
double ewd_i_short_energy(int i, const double *pos_array, const double *charge_array, int n_particles, double box_size)
{

    if (REAL_CUTOFF > box_size / 2)
    {
        LOG_WARNING("%s", CUTOFF_WARNING_MESSAGE);
    }

    double real_space_i_energy = 0;

    for (size_t j = 0; j < n_particles; j++)
    {

        // Exclude self particle in cell (0,0,0)
        if (i == j)
            continue;

        double r_ij_x = pos_array[c(i, 0)] - pos_array[c(j, 0)];
        double r_ij_y = pos_array[c(i, 1)] - pos_array[c(j, 1)];
        double r_ij_z = pos_array[c(i, 2)] - pos_array[c(j, 2)];

        // First image convention
        r_ij_x = pb_minimum_image(r_ij_x, box_size);
        r_ij_y = pb_minimum_image(r_ij_y, box_size);
        r_ij_z = pb_minimum_image(r_ij_z, box_size);

        const double r_ij_mod2 = r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z;

        // In first image convention _CUTOFF must be less than L/2
        if (r_ij_mod2 > REAL_CUTOFF * REAL_CUTOFF)
            continue;

        if (r_ij_mod2 < EWD_EPSILON)
        {
            printf("WARNING: 2 particles found below EWD_EPSILON distance, possible arithmetic error in energy evaluation");
        }

        const double r_ij_mod = sqrt(r_ij_mod2);

        real_space_i_energy += (charge_array[i] * charge_array[j]) * erfc(ALPHA * r_ij_mod) / r_ij_mod;
    }

    return real_space_i_energy;
}

/**
 * @brief Same behavior as "ewd_i_short_energy" but instead of looping over all the particles a VerletList is used.
 * The complexity is reduced from O(N) to O(1), where the constant depend on the mean number of neightbours stored in the verlet list
 * for the i particle.
 *
 * @note Requires a right initialized VerletList vl.
 * @note Complexity O(1)
 */
double ewd_verlet_i_short_energy(int i, const double *pos_array, const double *charge_array, const VerletList_t *vl, int n_particles, double box_size)
{

    if (REAL_CUTOFF > box_size / 2)
    {
        LOG_WARNING("%s", CUTOFF_WARNING_MESSAGE);
    }

    double real_space_i_energy = 0;

    for (size_t v = 0; v < vl[i].count; v++)
    {
        int j = vl[i].list[v];

        // Exclude self particle in cell (0,0,0)
        if (i == j)
            continue;

        double r_ij_x = pos_array[c(i, 0)] - pos_array[c(j, 0)];
        double r_ij_y = pos_array[c(i, 1)] - pos_array[c(j, 1)];
        double r_ij_z = pos_array[c(i, 2)] - pos_array[c(j, 2)];

        // First image convention
        r_ij_x = pb_minimum_image(r_ij_x, box_size);
        r_ij_y = pb_minimum_image(r_ij_y, box_size);
        r_ij_z = pb_minimum_image(r_ij_z, box_size);

        double r_ij_mod2 = r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z;

        // In first image convention _CUTOFF must be less than L/2
        if (r_ij_mod2 > REAL_CUTOFF * REAL_CUTOFF)
            continue;

        if (r_ij_mod2 < EWD_EPSILON)
        {
            printf("WARNING: 2 particles found below EWD_EPSILON distance, possible arithmetic error in energy evaluation");
        }

        double r_ij_mod = sqrt(r_ij_mod2);

        real_space_i_energy += (charge_array[i] * charge_array[j]) * erfc(ALPHA * r_ij_mod) / r_ij_mod;
    }

    return real_space_i_energy;
}

/**
 * @brief Compute the total real-space Coulomb looping over "ewd_i_short_energy". See that for more informations.
 *
 * @note Complexity of N*O(N) = O(N^2)
 */
double ewd_short_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    double real_space_energy = 0;

    for (size_t i = 0; i < n_particles; i++)
    {
        real_space_energy += ewd_i_short_energy(i, pos_array, charge_array, n_particles, box_size);
    }

    // Divide by 2 because in this implementation I am double counting
    return real_space_energy * 0.5;
}

/**
 * @brief Compute the total real-space Coulomb looping over "ewd_verlet_i_short_energy". See that for more informations.
 *
 * Complexity of N*O(1) = O(N)
 *
 * @note requires a right initialized VerletList vl.
 */
double ewd_verlet_short_energy(const double *pos_array, const double *charge_array, const VerletList_t *vl, int n_particles, double box_size)
{
    double real_space_energy = 0;

    for (size_t i = 0; i < n_particles; i++)
    {
        real_space_energy += ewd_verlet_i_short_energy(i, pos_array, charge_array, vl, n_particles, box_size);
    }

    // Divide by 2 because in this implementation I am double counting
    return real_space_energy * 0.5;
}

/**
 * @brief Compute the structural factor in reciprocal space for a given wavevector.
 *
 * This function calculates the structural factor S(k), which is the Fourier transform
 * of the charge distribution in the system, for a given reciprocal lattice vector (k = (k_x, k_y, k_z)).
 * The structural factor is used to calculate the long-range Coulomb interactions in reciprocal space in the Ewald summation.
 */
static inline double complex ewd_compute_structural_factor(const double *pos_array, const double *charge_array, int n_particles, double k_x, double k_y, double k_z)
{
    double complex structural_factor = 0;
    // Term O(N)
    for (size_t i = 0; i < n_particles; i++)
    {

        double r_i_x = pos_array[c(i, 0)];
        double r_i_y = pos_array[c(i, 1)];
        double r_i_z = pos_array[c(i, 2)];

        double dot_prod = (k_x * r_i_x + k_y * r_i_y + k_z * r_i_z);
        structural_factor += charge_array[i] * cexp(I * dot_prod);
    }

    return structural_factor;
}

/**
 * @brief Compute the reciprocal-space Coulomb energy of a system of charged particles.
 *
 * Calculates the long-range Coulomb interactions between particles in reciprocal space
 * (Fourier space) using the Ewald summation method. The interactions are summed over
 * reciprocal lattice vectors up to a specified range, with a screening factor determined
 * by the parameter `alpha`. This term captures the long-range interactions that decay
 * exponentially in reciprocal space.
 *
 * The Coulomb energy is calculated by summing over the structure factors for all
 * reciprocal lattice vectors within the defined range, applying the appropriate
 * screening factor, and normalizing by the system volume.
 */
double ewd_long_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    const double base_frequency = (2 * PI / box_size);
    double sum = 0;
    for (int k_x = -RECIPROCAL_RANGE; k_x <= RECIPROCAL_RANGE; k_x++)
    {
        for (int k_y = -RECIPROCAL_RANGE; k_y <= RECIPROCAL_RANGE; k_y++)
        {
            for (int k_z = -RECIPROCAL_RANGE; k_z <= RECIPROCAL_RANGE; k_z++)
            {
                // Ignore cell (0,0,0) in k-space
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;

                // K sphere
                if (k_x * k_x + k_y * k_y + k_z * k_z > RECIPROCAL_RANGE * RECIPROCAL_RANGE)
                    continue;

                double k_x_f = k_x * base_frequency;
                double k_y_f = k_y * base_frequency;
                double k_z_f = k_z * base_frequency;

                double k_mod2 = k_x_f * k_x_f + k_y_f * k_y_f + k_z_f * k_z_f;
                double complex structural_factor = ewd_compute_structural_factor(pos_array, charge_array, n_particles, k_x_f, k_y_f, k_z_f);
                sum += (structural_factor * conj(structural_factor)) * exp(-k_mod2 / (4 * ALPHA * ALPHA)) / k_mod2;
            }
        }
    }
    double volume = box_size * box_size * box_size;
    return 0.5 * (4 * PI) / volume * sum;
}

/**
 * @brief Compute the total Coulomb energy of a system of charged particles using the Ewald summation method.
 *
 * This function calculates the total Coulomb energy by summing the short-range real-space interactions,
 * the long-range reciprocal-space interactions, and the self-energy of the particles. The total energy
 * is computed by combining these components, subtracting the self-energy to avoid double-counting.
 *
 * The short-range interactions are computed using a cutoff in real space, and the long-range interactions
 * are computed using Fourier space. The self-energy accounts for the interaction of each particle with its
 * own charge distribution.
 *
 * @note Complexity O(N^2)
 */
double ewd_total_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{

    if (OPTIMIZED == 0)
    {
        printf("WARNING: Ewald parameters has not been optimized. Call \"optimizeParameter\" or set ALPHA, REAL_CUTOFF and RECIPROCAL_RANGE by hand.\n");
    }

    static double self_energy = 0;

    // O(N^2)
    double short_range_energy = ewd_short_energy(pos_array, charge_array, n_particles, box_size);
    // O(N)
    double long_range_energy = ewd_long_energy(pos_array, charge_array, n_particles, box_size);

    // O(1)
    // If the charge of the particle does not change the self_energy is constant so it has to be computed only once,
    // unless all particles are charge neutral self energy is always greater than 0
    if (self_energy == 0)
    {
        self_energy = ewd_self_energy(charge_array, n_particles);
    }

    return (short_range_energy + long_range_energy - self_energy);
}

/**
 * @brief Compute the total Coulomb energy like in "ewd_total_energy" but using a VerletList to reduce complexity.
 * Look at "ewd_total_energy" description for more informations.
 *
 * @note Complexity O(N)
 */
double ewd_verlet_total_energy(const double *pos_array, const double *charge_array, const VerletList_t *vl, int n_particles, double box_size)
{

    if (OPTIMIZED == 0)
    {
        printf("WARNING: Ewald parameters has not been optimized. Call \"optimizeParameter\" or set ALPHA, REAL_CUTOFF and RECIPROCAL_RANGE by hand.\n");
    }

    static double self_energy = 0;

    // O(N)
    double short_range_energy = ewd_verlet_short_energy(pos_array, charge_array, vl, n_particles, box_size);
    // O(N)
    double long_range_energy = ewd_long_energy(pos_array, charge_array, n_particles, box_size);

    // O(1)
    // If the charge of the particle does not change the self_energy is constant so it has to be computed only once,
    // unless all particles are charge neutral self energy is always greater than 0
    if (self_energy == 0)
    {
        self_energy = ewd_self_energy(charge_array, n_particles);
    }

    return (short_range_energy + long_range_energy - self_energy);
}

/**
 * @brief Compute the difference in the long energy given old_pos_array and new_pos_array where ONLY the i particle has moved
 * from r_i -> r_i'.
 *
 * The result of this function is exaclty the same as computing
 *
 * dE_long = ewd_long_energy(new_pos_array,...) - ewd_long_energy(old_pos_array,..)
 *
 * but it is faster since is does not compute the structural factor S(k), so it misses the loop over N bringing
 * the complexity from O(N) to O(1). S(K) is cmputed only once at the start of the simulation using the start postion array
 * and the function "ewd_init_S_K", then S_K has to be updated by hand using "ewd_update_S_K".
 * In the MC case we will update S(k) only if the step is accepted.
 * The update for S_K is also O(1).
 *
 * @note Complexity O(1)
 */
double ewd_delta_long_energy(int i, const double *new_pos_array, const double *old_pos_array, const double *charge_array, double box_size)
{

    const int NK = 2 * RECIPROCAL_RANGE + 1;
    const double base_frequency = (2 * PI / box_size);

    double delta_E = 0;

    for (int k_x = -RECIPROCAL_RANGE; k_x <= RECIPROCAL_RANGE; k_x++)
    {
        for (int k_y = -RECIPROCAL_RANGE; k_y <= RECIPROCAL_RANGE; k_y++)
        {
            for (int k_z = -RECIPROCAL_RANGE; k_z <= RECIPROCAL_RANGE; k_z++)
            {
                // Ignore cell (0,0,0) in k-space
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;

                // K sphere
                if (k_x * k_x + k_y * k_y + k_z * k_z > RECIPROCAL_RANGE * RECIPROCAL_RANGE)
                    continue;

                const double k_x_f = k_x * base_frequency;
                const double k_y_f = k_y * base_frequency;
                const double k_z_f = k_z * base_frequency;

                const double phase_old = k_x_f * old_pos_array[c(i, 0)] + k_y_f * old_pos_array[c(i, 1)] + k_z_f * old_pos_array[c(i, 2)];
                const double phase_new = k_x_f * new_pos_array[c(i, 0)] + k_y_f * new_pos_array[c(i, 1)] + k_z_f * new_pos_array[c(i, 2)];

                const int idx = (k_x + RECIPROCAL_RANGE) * NK * NK + (k_y + RECIPROCAL_RANGE) * NK + (k_z + RECIPROCAL_RANGE);

                const double complex dS_i = charge_array[i] * (cexp(I * phase_new) - cexp(I * phase_old));

                const double k2 = k_x_f * k_x_f + k_y_f * k_y_f + k_z_f * k_z_f;

                const double prefactor_k = (2 * PI / (box_size * box_size * box_size)) * exp(-k2 / (4 * ALPHA * ALPHA)) / k2;

                delta_E += prefactor_k * (2 * creal(conj(S_K[idx]) * dS_i) + creal(dS_i * conj(dS_i)));
            }
        }
    }
    return delta_E;
}

/**
 * @brief Update S_K given old_pos_array and new_pos_array where ONLY the i particle has moved
 * from r_i -> r_i'.
 * For more information look at the "ewd_delta_long_energy" description.
 *
 * @note Complexity O(1)
 */
void ewd_update_S_K(int i, const double *new_pos_array, const double *old_pos_array, const double *charge_array, double box_size)
{
    const int NK = 2 * RECIPROCAL_RANGE + 1;
    const double base_frequency = (2 * PI / box_size);

    for (int k_x = -RECIPROCAL_RANGE; k_x <= RECIPROCAL_RANGE; k_x++)
    {
        for (int k_y = -RECIPROCAL_RANGE; k_y <= RECIPROCAL_RANGE; k_y++)
        {
            for (int k_z = -RECIPROCAL_RANGE; k_z <= RECIPROCAL_RANGE; k_z++)
            {
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;
                if (k_x * k_x + k_y * k_y + k_z * k_z > RECIPROCAL_RANGE * RECIPROCAL_RANGE)
                    continue;

                const double k_x_f = k_x * base_frequency;
                const double k_y_f = k_y * base_frequency;
                const double k_z_f = k_z * base_frequency;

                const double phase_old = k_x_f * old_pos_array[c(i, 0)] + k_y_f * old_pos_array[c(i, 1)] + k_z_f * old_pos_array[c(i, 2)];
                const double phase_new = k_x_f * new_pos_array[c(i, 0)] + k_y_f * new_pos_array[c(i, 1)] + k_z_f * new_pos_array[c(i, 2)];

                const int idx = (k_x + RECIPROCAL_RANGE) * NK * NK + (k_y + RECIPROCAL_RANGE) * NK + (k_z + RECIPROCAL_RANGE);

                const double complex dS = charge_array[i] * (cexp(I * phase_new) - cexp(I * phase_old));

                S_K[idx] += dS;
            }
        }
    }
}

/**
 * @brief Used in "ewd_init_S_K", look at it's description for more inormation.
 */
static void ewd_fill_S_K(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    const int NK = 2 * RECIPROCAL_RANGE + 1;
    const double base_frequency = (2 * PI / box_size);

    for (int k_x = -RECIPROCAL_RANGE; k_x <= RECIPROCAL_RANGE; k_x++)
    {
        for (int k_y = -RECIPROCAL_RANGE; k_y <= RECIPROCAL_RANGE; k_y++)
        {
            for (int k_z = -RECIPROCAL_RANGE; k_z <= RECIPROCAL_RANGE; k_z++)
            {

                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;

                // K sphere of radius RECIPROCAL_RANGE in reciprocal lattice space
                if (k_x * k_x + k_y * k_y + k_z * k_z > RECIPROCAL_RANGE * RECIPROCAL_RANGE)
                    continue;

                const double k_x_f = k_x * base_frequency;
                const double k_y_f = k_y * base_frequency;
                const double k_z_f = k_z * base_frequency;

                const int idx = (k_x + RECIPROCAL_RANGE) * NK * NK + (k_y + RECIPROCAL_RANGE) * NK + (k_z + RECIPROCAL_RANGE);

                double complex Sk = 0.0 + 0.0 * I;

                for (int i = 0; i < n_particles; i++)
                {
                    const double phase =
                        k_x_f * pos_array[c(i, 0)] +
                        k_y_f * pos_array[c(i, 1)] +
                        k_z_f * pos_array[c(i, 2)];

                    Sk += charge_array[i] * cexp(I * phase);
                }

                S_K[idx] = Sk;
            }
        }
    }
}

/**
 * @brief Precompile and store in memory every S(K).
 * K_x assumes values in {-RECIPROCAL_RANGE, -RECIPROCAL_RANGE + 1,..., RECIPROCAL_RANGE - 1, RECIPROCAL_RANGE},
 * same for K_and K_z. Gived that, S(K) is stored as a flattered 3D matrix of size NK^3, where NK = 2 * RECIPROCAL_RANGE + 1;
 *
 * To retrieve the correct index for a vector K use:
 * idx = (k_x + RECIPROCAL_RANGE) * NK * NK + (k_y + RECIPROCAL_RANGE) * NK + (k_z + RECIPROCAL_RANGE);
 */
void ewd_init_S_K(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    const int NK = 2 * RECIPROCAL_RANGE + 1;
    const int Ntot = NK * NK * NK;

    if (S_K != NULL)
    {
        free(S_K);
    }

    S_K = calloc(Ntot, sizeof(double complex));
    if (!S_K)
        exit(EXIT_FAILURE);

    ewd_fill_S_K(pos_array, charge_array, n_particles, box_size);
}

#endif