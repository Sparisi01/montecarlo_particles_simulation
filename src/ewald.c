#ifndef EWALD_H
#define EWALD_H

#include <math.h>
#include <inttypes.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>

#include "constants.c"
#include "periodic_boundaries.c"

double _CUTOFF = (BOX_SIZE / 2);

/*
double errorsDifference(double error, double s, double Q, double cell_length, double alpha)
{
    return exp(-(s * s)) / (pow(s, 3. / 2)) * Q * sqrt((2 + PI) / (2 * PI * alpha * pow(cell_length, 3))) - error;
};

// Time constants
#define TAU_S 3.6
#define TAU_L 1.0
#define TAU_RAPP (TAU_S / TAU_L)

double findSbybisection(double a, double b, double error, double Q, double cell_length, double alpha, double precision)
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

    printf("ROOT FOUND: %.5E\n", c);
    return c;
}

double optimizeParameter(double error, double cell_length, int N_particles, double Q)
{
    double alpha = pow((TAU_RAPP * pow(PI, 3) / pow(cell_length, 6) * N_particles), 1. / 6);

    printf("ALPHA: %.5E\n", alpha);

    double s = findSbybisection(1e-3, 1e3, error, Q, cell_length, alpha, 1e-10);

    double rc = s / alpha;
    double kc = 2 * s * alpha;

    _ALPHA = alpha;
    _CUTOFF = rc;
    _K_RANGE_EWALD = ceil(kc / (2 * PI / cell_length));

    printf("OPTIMIZED PARAMETERS: R_C = %.5E, N_C = %.5E\n", rc, kc / (2 * PI / cell_length));
}
*/

// Parameter oprimization
static inline double ewd_alpha_by_precision(double precision)
{
    return sqrt(-log(precision));
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
static double ewd_self_energy(const double *charge_array, int n_particles, double alpha)
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
static double ewd_i_real_space_coulomb_energy(int i, double *pos_array, const double *charge_array, int n_particles, double box_size, double alpha)
{
    const int _R_RANGE_EWALD = 0; // Set to 0 first image convention

    if (_CUTOFF >= BOX_SIZE / 2)
    {
        printf("Cutoff to big for first image convention");
    }

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
                    if (r_ij_mod2 > _CUTOFF * _CUTOFF)
                        continue;
                    if (r_ij_mod2 == 0)
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

double ewd_real_space_coulomb_energy(double *pos_array, const double *charge_array, int n_particles, double box_size, double alpha)
{
    double real_space_energy = 0;

    for (size_t i = 0; i < n_particles; i++)
    {
        real_space_energy += ewd_i_real_space_coulomb_energy(i, pos_array, charge_array, n_particles, box_size, alpha);
    }

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
    for (size_t i = 0; i < n_particles; i++)
    {
        double qi = charge_array[i];

        double r_i_x = pos_array[c(i, 0)];
        double r_i_y = pos_array[c(i, 1)];
        double r_i_z = pos_array[c(i, 2)];

        double dot_prod = (k_x * r_i_x + k_y * r_i_y + k_z * r_i_z);
        structural_factor += qi * cexp(I * dot_prod);
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
static double ewd_reciprocal_space_coulomb_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size, double alpha)
{

    const double base_frequency = (2 * PI / box_size);
    const int _K_RANGE_EWALD = 8;

    double sum = 0;
    for (int k_x = -_K_RANGE_EWALD; k_x <= _K_RANGE_EWALD; k_x++)
    {
        for (int k_y = -_K_RANGE_EWALD; k_y <= _K_RANGE_EWALD; k_y++)
        {
            for (int k_z = -_K_RANGE_EWALD; k_z <= _K_RANGE_EWALD; k_z++)
            {
                // Ignore cell (0,0,0) in k-space
                if (k_x == 0 && k_y == 0 && k_z == 0)
                    continue;

                // K sphere
                if (k_x * k_x + k_y * k_y + k_z * k_z > _K_RANGE_EWALD * _K_RANGE_EWALD)
                    continue;

                double k_x = k_x * base_frequency;
                double k_y = k_y * base_frequency;
                double k_z = k_z * base_frequency;

                double k_mod2 = k_x * k_x + k_y * k_y + k_z * k_z;
                double complex structural_factor = edw_compute_structural_factor(pos_array, charge_array, n_particles, k_x, k_y, k_z);
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
double ewd_total_energy(const double *pos_array, const double *charge_array, int n_particles, double box_size)
{
    double alpha = ewd_alpha_by_precision(1e-6);

    double short_range_energy = ewd_real_space_coulomb_energy(pos_array, charge_array, n_particles, box_size, alpha);
    double long_range_energy = ewd_reciprocal_space_coulomb_energy(pos_array, charge_array, n_particles, box_size, alpha);
    double self_energy = ewd_self_energy(charge_array, n_particles, alpha);

    return (short_range_energy + long_range_energy - self_energy);
}

#endif