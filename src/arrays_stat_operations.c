#ifndef ARRAYS_STAT_OPERATIONS_C
#define ARRAYS_STAT_OPERATIONS_C

#include <stdlib.h>
#include "logger.c"

int pointerValidator(const void *ptr)
{
    if (ptr == NULL)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

/**
 * @brief given an array arr of double and size N compute the mean = sum(arr[i])/N
 */
double array_mean(const double *array, unsigned int array_dim)
{
    if (!pointerValidator(array))
    {
        LOG_FATAL("Array null pointer");
    }

    if (array_dim == 0)
    {
        LOG_WARNING("array_dim equal to 0, 0 returned");
        return 0;
    }

    double mean = 0;
    for (size_t i = 0; i < array_dim; i++)
    {
        mean += array[i];
    }
    return mean / array_dim;
}

/**
 * @brief given an array arr of double and size N compute the mean of squares = sum(arr[i]^2)/N
 */
double array_mean2(const double *array, unsigned int array_dim)
{
    if (!pointerValidator(array))
    {
        LOG_FATAL("Array null pointer");
    }

    if (array_dim == 0)
    {
        LOG_WARNING("array_dim equal to 0, 0 returned");
        return 0;
    }

    double mean2 = 0;
    for (size_t i = 0; i < array_dim; i++)
    {
        mean2 += array[i] * array[i];
    }
    return mean2 / array_dim;
}

/**
 * @brief given an array arr of double and size N compute the variance
 */
double array_var(const double *array, unsigned int array_dim)
{
    if (!pointerValidator(array))
    {
        LOG_FATAL("Array null pointer");
    }

    if (array_dim == 0)
    {
        LOG_WARNING("array_dim equal to 0, 0 returned");
        return 0;
    }

    return array_mean2(array, array_dim) - array_mean(array, array_dim) * array_mean(array, array_dim);
}

/**
 * @brief multiply each element of arr by the constant c
 */
void array_const_mult(double *array, unsigned int array_dim, double c)
{
    if (!pointerValidator(array))
    {
        LOG_FATAL("Array null pointer");
    }

    if (array_dim == 0)
    {
        LOG_WARNING("array_dim equal to 0, 0 returned");
        return;
    }

    for (size_t i = 0; i < array_dim; i++)
    {
        array[i] *= c;
    }
}

/**
 * @brief multiply each element of arr by the constant c
 */
double array_dot_product(const double *array_1, const double *array_2, unsigned int array_dim)
{
    if (!pointerValidator(array_1))
    {
        LOG_FATAL("Array_1 null pointer");
    }

    if (!pointerValidator(array_2))
    {
        LOG_FATAL("Array_2 null pointer");
    }

    if (array_dim == 0)
    {
        LOG_WARNING("array_dim equal to 0, 0 returned");
        return 0;
    }

    double sum = 0;
    for (size_t i = 0; i < array_dim; i++)
    {
        sum += array_1[i] * array_2[i];
    }
    return sum;
}

#endif