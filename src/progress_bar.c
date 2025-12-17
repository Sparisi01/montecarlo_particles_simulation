#ifndef PROGRESS_BAR
#define PROGRESS_BAR

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

/* PRINT_PROGRESS
 *
 * Prints a terminal progress bar with percentage and estimated time remaining (ETA).
 *
 * Parameters:
 *   current     – current iteration (0 … total)
 *   total       – total number of iterations
 *   start_time  – the time (time_t) when the process started; used to compute ETA
 *
 * Behavior:
 *   - The function prints a character progress bar using "#" character.
 *   - It overwrites the same terminal line using '\r'.
 *   - ETA is computed from elapsed time and current progress and displayed as HH:MM:SS.
 *
 * Usage:
 *   Call this function periodically inside a loop to update the progress bar without
 *   significantly impacting performance.
 */
void print_progress(size_t current, size_t total, clock_t start_time)
{
    const int barWidth = 20;
    double progress = (double)current / total;
    int filled = progress * barWidth;

    // Time elapsed
    clock_t now = clock();
    double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;

    // ETA stimation
    double eta = (progress > 0.0) ? elapsed * (1.0 - progress) / progress : 0.0;

    // Conversion of ETA hh:mm:ss
    int eta_h = (int)eta / 3600;
    int eta_m = ((int)eta % 3600) / 60;
    int eta_s = (int)eta % 60;
    printf("\r\033[2K[");

    for (int i = 0; i < barWidth; i++)
    {
        if (i <= filled)
            printf("#");
        else
            printf(".");
    }

    printf("] %5.1f%%  ETA: %02d:%02d:%02d (hh:mm:ss)",
           progress * 100, eta_h, eta_m, eta_s);

    fflush(stdout);
}

#endif