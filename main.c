#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/structures.c"
#include "src/constants.c"

int main(int argc, char const *argv[])
{
    // Set seed for reproducibility
    srand(SEED);

    struct particle *particle_array = (struct particle *)malloc(N * sizeof(struct particle));
    if (particle_array == NULL)
        exit(EXIT_FAILURE);
    
    // Init positions using uniform distribution in the square_box
    for (size_t i = 0; i < N; i++)
    {
        particle_array[i].pos.x = rand() / (RAND_MAX - 1.) * BOX_SIZE;
        particle_array[i].pos.y = rand() / (RAND_MAX - 1.) * BOX_SIZE;
        particle_array[i].pos.z = rand() / (RAND_MAX - 1.) * BOX_SIZE;
    }
    
    clock_t begin = clock();

    /*
        EVERYTHING WE WILL WRITE 
    */
    
    clock_t end = clock();

    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;

    printf("Total time: %.0lf ms", time_spent * 1000);


    free(particle_array);

    return 0;
}
