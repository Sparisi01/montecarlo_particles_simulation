#ifndef PERIODIC_BOUNDARIES_C
#define PERIODIC_BOUNDARIES_C

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/**
 * @brief Apply periodic boundary conditions to a position.
 *
 * @param x  Position coordinate to be wrapped
 * @param box_size
 * @return   Wrapped position in the range [0, L)
 */
// ANCHOR - succesfully tested with 1 particle and constant force along x axis
static inline double pb_wrap_position(double x, double box_size)
{
    return x - floor(x / box_size) * box_size;
}

/**
 * @brief Apply the minimum image convention on position component dx.

 * @param dx Displacement along x
 * @param box_size
 * @return   Minimum-image
 */
// ANCHOR - succesfully tested with 2 particle displaced in the corner of the box
static inline double pb_minimum_image(double dx, double box_size)
{
    dx -= floor(dx / box_size + 0.5) * box_size;
    return dx;
}

#endif