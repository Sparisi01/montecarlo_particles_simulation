#ifndef ARRAYS_STAT_OPERATIONS_C
#define ARRAYS_STAT_OPERATIONS_C

#include <stdlib.h>

/**
 * @brief given an array arr of double and size N compute the mean = sum(arr[i])/N
 */
double array_mean(double *arr, int n)
{
    double mean = 0;
    for (size_t i = 0; i < n; i++)
    {
        mean += arr[i];
    }
    return mean / n;
}

/**
 * @brief given an array arr of double and size N compute the mean of squares = sum(arr[i]^2)/N
 */
double array_mean2(double *arr, int n)
{
    double mean2 = 0;
    for (size_t i = 0; i < n; i++)
    {
        mean2 += arr[i] * arr[i];
    }
    return mean2 / n;
}

/**
 * @brief given an array arr of double and size N compute the variance
 */
double array_var(double *arr, int n)
{
    return array_mean2(arr, n) - array_mean(arr, n) * array_mean(arr, n);
}

#endif