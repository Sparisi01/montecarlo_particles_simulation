#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define N 100000
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

// CUDA kernel to compute energy
__global__ void gpu_compute_energy(particle* p_array, float* energy, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Only proceed if the thread index is within the bounds of the particle array
    if (idx < n) {
        float local_energy = 0.0;

        for (int j = idx + 1; j < n; j++) {
            float dx = p_array[idx].pos.x - p_array[j].pos.x;
            float dy = p_array[idx].pos.y - p_array[j].pos.y;
            float dz = p_array[idx].pos.z - p_array[j].pos.z;

            // Compute squared distance
            float distance_sq = dx * dx + dy * dy + dz * dz;
            local_energy += 1.0 / distance_sq;
        }

        // Use atomicAdd to safely accumulate the energy on the device
        atomicAdd(energy, local_energy);
    }
}

// Host function to compute energy using CUDA
float compute_energy_cuda(particle* p_array, int n) {
    particle* d_p_array;
    float* d_energy;
    float h_energy = 0.0;

    // Allocate memory on the device
    cudaMalloc((void**)&d_p_array, n * sizeof(particle));
    cudaMalloc((void**)&d_energy, sizeof(float));

    // Initialize the energy to zero on the device
    cudaMemset(d_energy, 0, sizeof(float));

    // Copy particles array from host to device
    cudaMemcpy(d_p_array, p_array, n * sizeof(particle), cudaMemcpyHostToDevice);

    // Launch kernel with enough blocks and threads
    int blockSize = 512;
    int numBlocks = (n + blockSize - 1) / blockSize;
    gpu_compute_energy<<<numBlocks, blockSize>>>(d_p_array, d_energy, n);

    // Check for errors in kernel execution
    cudaDeviceSynchronize();

    // Copy result from device to host
    cudaMemcpy(&h_energy, d_energy, sizeof(float), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_p_array);
    cudaFree(d_energy);

    return h_energy;
}

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

    float energy = compute_energy_cuda(particle_array, N);
    //float energy = compute_energy(particle_array, N);
    clock_t end = clock();
    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;


    printf("Energy: %lf\n",energy);
    printf("Time: %.0lf ms", time_spent*1000);
    free(particle_array);
    
    return 0;
}


//6992.375000
//6758.898926
