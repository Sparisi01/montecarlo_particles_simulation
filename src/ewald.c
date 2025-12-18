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

#define EWD_EPSILON 1e-8
#define EWD_PRECISION 1e-1

static double reciprocal_range;
static double real_cutoff;
static double alpha;
static int optimized = 0;

static double errorsDifference(double error, double s, double Q, double cell_length, double alpha)
{
    return exp(-(s * s)) / (pow(s, 3. / 2)) * Q * sqrt((2 + PI) / (2 * PI * alpha * pow(cell_length, 3))) - error;
};

// Time constants
#define TAU_S 3.6
#define TAU_L 1.0
#define TAU_RAPP (TAU_S / TAU_L)

static double findSbybisection(double a, double b, double error, double Q, double cell_length, double alpha, double precision)
{
    double max = errorsDifference(error, a, Q, cell_length, alpha);
    double min = errorsDifference(error, b, Q, cell_length, alpha);

    if (!(max > 0 && min < 0))
    {
        printf("Max e min non rispettano parametri bisezione\n");
    };
    double c = 0;
    int root_find = 0;
    while (!root_find)
    {
        c = (a + b) / 2;
        double mid_value = errorsDifference(error, c, Q, cell_length, alpha);
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

    // printf("ROOT FOUND: %.5E\n", c);
    return c;
}

void optimizeParameter(double error, double box_size, const double *charge_array, int n_particles)
{
    optimized = 1;

    double Q2 = 0;
    for (size_t i = 0; i < n_particles; i++)
    {
        Q2 += charge_array[i] * charge_array[i];
    }

    alpha = pow((TAU_RAPP * pow(PI, 3) / pow(box_size, 6) * n_particles), 1. / 6);

    double s = findSbybisection(1e-3, 1e3, error, Q2, box_size, alpha, 1e-10);
    real_cutoff = s / alpha;

    double kc = 2 * s * alpha;
    reciprocal_range = ceil(kc / (2 * PI / box_size));

    printf("OPTIMIZED PARAMETERS: ALPHA = %.5E, R_C = %.5E, N_C_K = %.5E\n", alpha, real_cutoff, reciprocal_range);

    if (real_cutoff > box_size / 2)
    {
        printf("WARNING: Cutoff to big for first image convention\n");
    }
}

// Simple parameter oprimization
static inline void ewd_alpha_by_precision(double precision)
{
    alpha = sqrt(-log(precision));
}

/**
 * @brief Compute the self-energy of a system of charged particles.
 *
 * Calculates the self-energy contribution of the Coulomb interaction for each particle
 * in the system. This term represents the energy of a particle due to its own charge
 * distribution in the system. The self-energy is computed by summing the squares of
 * the charges of all particles, applying the screening parameter `alpha`, and normalizing
 * by a factor of `1 / √π`.
 *
 * This function is typically used in the context of Ewald summation to handle the
 * self-interaction term. In the final energy calculation, this term is usually
 * excluded or adjusted to avoid double-counting.
 *
 * @param charge_array  Array of particle charges
 * @param n_particles   Total number of particles in the system
 * @param alpha         Screening parameter used in the Ewald summation
 *
 * @return The self-energy of the system, adjusted by the parameter `alpha` and normalized
 *         by `1 / √π`.
 *
 * @note The self-energy term is typically subtracted from the total energy to avoid
 *       double-counting in the pairwise interaction calculation.
 */
double ewd_self_energy(const double *charge_array, int n_particles)
{
    double sum = 0;
    for (size_t i = 0; i < n_particles; i++)
    {
        sum += charge_array[i] * charge_array[i];
    }
    return sum * alpha / SQRT_PI;
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
 * @param pos_array     Array of particle positions (flattened, 3 values per particle)
 * @param charge_array  Array of particle charges
 * @param n_particles   Total number of particles in the system
 * @param box_size      Size of the periodic simulation box (assumed cubic)
 * @param alpha         Screening parameter used in the real-space Coulomb potential
 *
 * @return The real-space Coulomb energy of the system, including interactions
 *         between particles within the cutoff distance, and screened by `alpha`.
 *
 * @warning If `_CUTOFF` parameter is invalid (<= 0), the function
 *          will terminate with an error message.
 */
double ewd_i_real_space_coulomb_energy(int i, const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    const int _R_RANGE_EWALD = 0; // Set to 0 first image convention

    double real_space_i_energy = 0;

    for (size_t j = 0; j < n_particles; j++)
    {
        for (int r_x = -_R_RANGE_EWALD; r_x <= _R_RANGE_EWALD; r_x++)
        {
            for (int r_y = -_R_RANGE_EWALD; r_y <= _R_RANGE_EWALD; r_y++)
            {
                for (int r_z = -_R_RANGE_EWALD; r_z <= _R_RANGE_EWALD; r_z++)
                {
                    // Exclude self particle in cell (0,0,0)
                    if (r_x == 0 && r_y == 0 && r_z == 0 && i == j)
                        continue;

                    double r_ij_x = pos_array[c(i, 0)] - pos_array[c(j, 0)];
                    double r_ij_y = pos_array[c(i, 1)] - pos_array[c(j, 1)];
                    double r_ij_z = pos_array[c(i, 2)] - pos_array[c(j, 2)];

                    // First image convention
                    r_ij_x = pb_minimum_image(r_ij_x, box_size) + r_x * box_size;
                    r_ij_y = pb_minimum_image(r_ij_y, box_size) + r_y * box_size;
                    r_ij_z = pb_minimum_image(r_ij_z, box_size) + r_z * box_size;

                    double r_ij_mod2 = r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z;

                    // In first image convention _CUTOFF must be less than L/2
                    if (r_ij_mod2 > real_cutoff * real_cutoff)
                        continue;
                    if (r_ij_mod2 < EWD_EPSILON)
                        goto PARTICLE_OVERLAP_ERROR;

                    // Avoid Sqrt if not needed
                    double r_ij_mod = sqrt(r_ij_mod2);

                    real_space_i_energy += (charge_array[i] * charge_array[j]) * erfc(alpha * r_ij_mod) / r_ij_mod;
                }
            }
        }
    }

    return real_space_i_energy;

// Error handling
PARTICLE_OVERLAP_ERROR:
    printf("ERROR: particle overlap");
    exit(EXIT_FAILURE);
}

double ewd_verlet_i_real_space_coulomb_energy(int i, const double *pos_array, const double *charge_array, const VerletList_t *vl, int n_particles, double box_size)
{
    const int _R_RANGE_EWALD = 0; // Set to 0 first image convention

    double real_space_i_energy = 0;

    for (size_t v = 0; v < vl[i].count; v++)
    {
        int j = vl[i].list[v];

        for (int r_x = -_R_RANGE_EWALD; r_x <= _R_RANGE_EWALD; r_x++)
        {
            for (int r_y = -_R_RANGE_EWALD; r_y <= _R_RANGE_EWALD; r_y++)
            {
                for (int r_z = -_R_RANGE_EWALD; r_z <= _R_RANGE_EWALD; r_z++)
                {
                    // Exclude self particle in cell (0,0,0)
                    if (r_x == 0 && r_y == 0 && r_z == 0 && i == j)
                        continue;

                    double r_ij_x = pos_array[c(i, 0)] - pos_array[c(j, 0)];
                    double r_ij_y = pos_array[c(i, 1)] - pos_array[c(j, 1)];
                    double r_ij_z = pos_array[c(i, 2)] - pos_array[c(j, 2)];

                    // First image convention
                    r_ij_x = pb_minimum_image(r_ij_x, box_size) + r_x * box_size;
                    r_ij_y = pb_minimum_image(r_ij_y, box_size) + r_y * box_size;
                    r_ij_z = pb_minimum_image(r_ij_z, box_size) + r_z * box_size;

                    double r_ij_mod2 = r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z;

                    // In first image convention _CUTOFF must be less than L/2
                    if (r_ij_mod2 > real_cutoff * real_cutoff)
                        continue;
                    if (r_ij_mod2 < EWD_EPSILON)
                        goto PARTICLE_OVERLAP_ERROR;

                    // Avoid Sqrt if not needed
                    double r_ij_mod = sqrt(r_ij_mod2);

                    real_space_i_energy += (charge_array[i] * charge_array[j]) * erfc(alpha * r_ij_mod) / r_ij_mod;
                }
            }
        }
    }

    return real_space_i_energy;

// Error handling
PARTICLE_OVERLAP_ERROR:
    printf("ERROR: particle overlap");
    exit(EXIT_FAILURE);
}

// NOTE - I preferred to separate the energy into i-energies for a possible future implementation of verlet list
//  in that case instead of passing pos_array one can just pass the verlet array associated with particle i
double ewd_real_space_coulomb_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    double real_space_energy = 0;

    for (size_t i = 0; i < n_particles; i++)
    {
        real_space_energy += ewd_i_real_space_coulomb_energy(i, pos_array, charge_array, n_particles, box_size);
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
 *
 * @param pos_array     Array of particle positions (flattened, 3 values per particle)
 * @param charge_array  Array of particle charges
 * @param n_particles   Total number of particles in the system
 * @param k_x           The x-component of the reciprocal lattice vector
 * @param k_y           The y-component of the reciprocal lattice vector
 * @param k_z           The z-component of the reciprocal lattice vector
 *
 * @return A complex number representing the structural factor S(k).
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
 *
 * @param pos_array     Array of particle positions (flattened, 3 values per particle)
 * @param charge_array  Array of particle charges
 * @param n_particles   Total number of particles in the system
 * @param box_size      Size of the periodic simulation box (assumed cubic)
 * @param alpha         Screening parameter used in the reciprocal-space Coulomb potential
 *
 * @return The reciprocal-space Coulomb energy of the system, including long-range
 *         interactions, screened by `alpha`.
 */
double ewd_reciprocal_space_coulomb_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    const double base_frequency = (2 * PI / box_size);
    double sum = 0;
    for (int k_x = -reciprocal_range; k_x <= reciprocal_range; k_x++)
    {
        for (int k_y = -reciprocal_range; k_y <= reciprocal_range; k_y++)
        {
            for (int k_z = -reciprocal_range; k_z <= reciprocal_range; k_z++)
            {
                // Ignore cell (0,0,0) in k-space
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;

                // K sphere
                if (k_x * k_x + k_y * k_y + k_z * k_z > reciprocal_range * reciprocal_range)
                    continue;

                double k_x_f = k_x * base_frequency;
                double k_y_f = k_y * base_frequency;
                double k_z_f = k_z * base_frequency;

                double k_mod2 = k_x_f * k_x_f + k_y_f * k_y_f + k_z_f * k_z_f;
                double complex structural_factor = ewd_compute_structural_factor(pos_array, charge_array, n_particles, k_x_f, k_y_f, k_z_f);
                sum += (structural_factor * conj(structural_factor)) * exp(-k_mod2 / (4 * alpha * alpha)) / k_mod2;
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
 * @param pos_array     Array of particle positions (flattened, 3 values per particle)
 * @param charge_array  Array of particle charges
 * @param n_particles   Total number of particles in the system
 * @param box_size      Size of the periodic simulation box (assumed cubic)
 *
 * @return The total Coulomb energy of the system, calculated as the sum of short-range interactions,
 *         long-range interactions, and self-energy (with the self-energy subtracted to avoid double-counting).
 *
 */
double ewd_total_coulomb_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    // double alpha = ewd_alpha_by_precision(EWD_PRECISION);

    if (optimized == 0)
    {
        optimizeParameter(EWD_PRECISION, box_size, charge_array, n_particles);
    }

    static double self_energy = 0;

    // O(N^2)
    double short_range_energy = ewd_real_space_coulomb_energy(pos_array, charge_array, n_particles, box_size);
    // O(N)
    double long_range_energy = ewd_reciprocal_space_coulomb_energy(pos_array, charge_array, n_particles, box_size);

    // O(1)
    // If the charge of the particle does not change the self_energy is constant so it has to be computed only once,
    // unless all particles are charge neutral self energy is always greater than 0
    if (self_energy == 0)
    {
        self_energy = ewd_self_energy(charge_array, n_particles);
    }

    return (short_range_energy + long_range_energy - self_energy);
}

static double complex *S_k = NULL;
static int Nk_global;

static void fillS_k(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    const double base_frequency = (2 * PI / box_size);

    for (int k_x = -reciprocal_range; k_x <= reciprocal_range; k_x++)
    {
        for (int k_y = -reciprocal_range; k_y <= reciprocal_range; k_y++)
        {
            for (int k_z = -reciprocal_range; k_z <= reciprocal_range; k_z++)
            {
                // Ignore cell (0,0,0) in k-space
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;

                // K sphere
                if (k_x * k_x + k_y * k_y + k_z * k_z > reciprocal_range * reciprocal_range)
                    continue;

                double k_x_f = k_x * base_frequency;
                double k_y_f = k_y * base_frequency;
                double k_z_f = k_z * base_frequency;

                int idx = (k_x + reciprocal_range) * Nk_global * Nk_global + (k_y + reciprocal_range) * Nk_global + (k_z + reciprocal_range);

                double complex Sk = 0.0 + 0.0 * I;

                for (int i = 0; i < n_particles; i++)
                {
                    double phase =
                        k_x_f * pos_array[c(i, 0)] +
                        k_y_f * pos_array[c(i, 1)] +
                        k_z_f * pos_array[c(i, 2)];

                    Sk += charge_array[i] * cexp(I * phase);
                }

                S_k[idx] = Sk;
            }
        }
    }
}

void init_Sk(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    Nk_global = 2 * reciprocal_range + 1;
    int Ntot = Nk_global * Nk_global * Nk_global;

    S_k = calloc(Ntot, sizeof(double complex));
    if (!S_k)
        exit(EXIT_FAILURE);

    fillS_k(pos_array, charge_array, n_particles, box_size);
}

double ewd_delta_reciprocal_energy(int i, const double *new_pos_array, const double *old_pos_array, const double *charge_array, double box_size)
{

    const double base_frequency = (2 * PI / box_size);

    double delta_E = 0;

    for (int k_x = -reciprocal_range; k_x <= reciprocal_range; k_x++)
    {
        for (int k_y = -reciprocal_range; k_y <= reciprocal_range; k_y++)
        {
            for (int k_z = -reciprocal_range; k_z <= reciprocal_range; k_z++)
            {
                // Ignore cell (0,0,0) in k-space
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;

                // K sphere
                if (k_x * k_x + k_y * k_y + k_z * k_z > reciprocal_range * reciprocal_range)
                    continue;

                double k_x_f = k_x * base_frequency;
                double k_y_f = k_y * base_frequency;
                double k_z_f = k_z * base_frequency;

                double phase_old = k_x_f * old_pos_array[c(i, 0)] + k_y_f * old_pos_array[c(i, 1)] + k_z_f * old_pos_array[c(i, 2)];
                double phase_new = k_x_f * new_pos_array[c(i, 0)] + k_y_f * new_pos_array[c(i, 1)] + k_z_f * new_pos_array[c(i, 2)];

                int idx = (k_x + reciprocal_range) * Nk_global * Nk_global + (k_y + reciprocal_range) * Nk_global + (k_z + reciprocal_range);

                double complex dS_i = charge_array[i] * (cexp(I * phase_new) - cexp(I * phase_old));

                double k2 = k_x_f * k_x_f + k_y_f * k_y_f + k_z_f * k_z_f;

                double prefactor_k = (2 * PI / (box_size * box_size * box_size)) * exp(-k2 / (4 * alpha * alpha)) / k2;

                delta_E += prefactor_k * (2 * creal(conj(S_k[idx]) * dS_i) + creal(dS_i * conj(dS_i)));
            }
        }
    }
    return delta_E;
}

void updateS_k(int i, const double *new_pos_array, const double *old_pos_array, const double *charge_array, double box_size)
{
    const double base_frequency = (2 * PI / box_size);
    for (int k_x = -reciprocal_range; k_x <= reciprocal_range; k_x++)
    {
        for (int k_y = -reciprocal_range; k_y <= reciprocal_range; k_y++)
        {
            for (int k_z = -reciprocal_range; k_z <= reciprocal_range; k_z++)
            {
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;
                if (k_x * k_x + k_y * k_y + k_z * k_z > reciprocal_range * reciprocal_range)
                    continue;

                double k_x_f = k_x * base_frequency;
                double k_y_f = k_y * base_frequency;
                double k_z_f = k_z * base_frequency;

                double phase_old = k_x_f * old_pos_array[c(i, 0)] + k_y_f * old_pos_array[c(i, 1)] + k_z_f * old_pos_array[c(i, 2)];
                double phase_new = k_x_f * new_pos_array[c(i, 0)] + k_y_f * new_pos_array[c(i, 1)] + k_z_f * new_pos_array[c(i, 2)];

                int idx = (k_x + reciprocal_range) * Nk_global * Nk_global + (k_y + reciprocal_range) * Nk_global + (k_z + reciprocal_range);

                double complex dS = charge_array[i] * (cexp(I * phase_new) - cexp(I * phase_old));
                S_k[idx] += dS;
            }
        }
    }
}

#endif