#ifndef PERIODIC_BOUNDARIES
#define PERIODIC_BOUNDARIES

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/**
 * @brief Apply periodic boundary conditions to a position.
 *
 * Wraps a coordinate `x` into the primary simulation box of size `L`
 * using periodic boundary conditions.
 *
 * The returned value is guaranteed to lie in the interval:
 *     0 ≤ x < L
 *
 * This function is intended for wrapping particle positions, not
 * displacement vectors. It correctly handles values of `x` that are
 * outside the box by an arbitrary number of box lengths, including
 * negative values.
 *
 * Mathematical form:
 *     x_wrapped = x - floor(x / L) * L
 *
 * @param x  Position coordinate to be wrapped
 * @param box_size  Size of the periodic box (must be > 0)
 * @return   Wrapped position in the range [0, L)
 */
// ANCHOR - succesfully tested with 1 particle and constant force along x axis
static inline double pb_wrap_position(double x, double box_size)
{
    return x - floor(x / box_size) * box_size;
}

/**
 * @brief Apply the minimum image convention on position component dx.
 *
 * Given the displacement position dx between two particles,
 * applies periodic boundary conditions so that dx lies
 * in the range [-L/2, L/2].
 *
 * @param dx Displacement along x
 * @param box_size  Simulation box size (assumed cubic, L > 0)
 * @return   Minimum-image
 */
// ANCHOR - succesfully tested with 2 particle displaced in the corner of the box
static inline double pb_minimum_image(double dx, double box_size)
{
    dx -= floor(dx / box_size + 0.5) * box_size;
    return dx;
}

#endif