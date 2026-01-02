#ifndef LOGGER_C
#define LOGGER_C

#include <stdio.h>
#include <stdlib.h>

#define COLOR_RED "\033[31m"
#define COLOR_YELLOW "\033[33m"
#define COLOR_BLUE "\033[34m"
#define COLOR_GREEN "\033[32m"
#define COLOR_RESET "\033[0m"

#define STYLE_BOLD "\033[1m"
#define STYLE_RESET "\033[0m"

#define LOG_ERROR(fmt, ...)                                                                  \
    fprintf(stderr, STYLE_BOLD COLOR_RED "ERROR (%s:%d): " COLOR_RESET STYLE_RESET fmt "\n", \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define LOG_WARNING(fmt, ...)                                                                     \
    fprintf(stderr, STYLE_BOLD COLOR_YELLOW "WARNING (%s:%d): " COLOR_RESET STYLE_RESET fmt "\n", \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define LOG_INFO(fmt, ...)                                                           \
    fprintf(stderr, STYLE_BOLD COLOR_BLUE "INFO: " STYLE_RESET COLOR_RESET fmt "\n", \
            ##__VA_ARGS__)

#define LOG_TEST(fmt, ...)                                                            \
    fprintf(stderr, STYLE_BOLD COLOR_GREEN "TEST: " STYLE_RESET COLOR_RESET fmt "\n", \
            ##__VA_ARGS__)

#define LOG_FATAL(fmt, ...)            \
    do                                 \
    {                                  \
        LOG_ERROR(fmt, ##__VA_ARGS__); \
        exit(EXIT_FAILURE);            \
    } while (0)

#endif
