#ifndef LOGGER_C
#define LOGGER_C

#include <stdio.h>
#include <stdlib.h>

#ifdef LOGGER_NO_COLOR
#define COLOR_RED ""
#define COLOR_YELLOW ""
#define COLOR_BLUE ""
#define COLOR_RESET ""
#else
#define COLOR_RED "\033[31m"
#define COLOR_YELLOW "\033[33m"
#define COLOR_BLUE "\033[34m"
#define COLOR_RESET "\033[0m"
#endif

#define LOG_ERROR(fmt, ...)                                           \
    fprintf(stderr, COLOR_RED "ERROR (%s:%d): " fmt COLOR_RESET "\n", \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define LOG_WARNING(fmt, ...)                                              \
    fprintf(stderr, COLOR_YELLOW "WARNING (%s:%d): " fmt COLOR_RESET "\n", \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define LOG_INFO(fmt, ...)             \
    fprintf(stderr, "INFO: " fmt "\n", \
            ##__VA_ARGS__)

#define LOG_FATAL(fmt, ...)            \
    do                                 \
    {                                  \
        LOG_ERROR(fmt, ##__VA_ARGS__); \
        exit(EXIT_FAILURE);            \
    } while (0)

#endif
