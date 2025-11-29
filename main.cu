#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define N 100
#define SEED 42
#define BOX_SIZE 100

struct vec3 {
    float x,y,z;
} typedef vec3;

struct particle {
    vec3 pos;
    vec3 vel;
    int charge;
    float mass;
} typedef particle;

float compute_energy(particle* p_array, int n){
    float energy = 0;
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i+1; j < n; j++)
        {
           float dx_sq = pow(p_array[i].pos.x-p_array[j].pos.x,2);
           float dy_sq = pow(p_array[i].pos.y-p_array[j].pos.y,2);
           float dz_sq = pow(p_array[i].pos.z-p_array[j].pos.z,2);
           energy += 1/(dx_sq+dy_sq+dz_sq);
        }
    }

    return energy;
}

int main(int argc, char const *argv[])
{
    // Set seed for reproducibility
    srand(SEED);

    particle* particle_array = (particle*)malloc(N*sizeof(particle));

    for (size_t i = 0; i < N; i++)
    {
        particle_array[i].pos.x = rand()/(RAND_MAX-1.l) * BOX_SIZE;
        particle_array[i].pos.y = rand()/(RAND_MAX-1.l) * BOX_SIZE;
        particle_array[i].pos.z = rand()/(RAND_MAX-1.l) * BOX_SIZE;
    }
    
    clock_t begin = clock();

    float energy = compute_energy(particle_array, N);

    clock_t end = clock();

    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;


    printf("Energy: %lf\n",energy);
    printf("Time: %.0lf ms", time_spent*1000);
    free(particle_array);
    
    return 0;
}


//6992.375000
//6758.898926
