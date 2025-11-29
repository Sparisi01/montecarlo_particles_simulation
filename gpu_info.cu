#include <stdio.h>
#include <cuda_runtime.h>

int main() {
    int device_count;
    cudaGetDeviceCount(&device_count);
    printf("Numero di GPU disponibili: %d\n", device_count);

    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp device_properties;
        cudaGetDeviceProperties(&device_properties, i);

        printf("\nProprietà della GPU %d:\n", i);
        printf("Nome: %s\n", device_properties.name);
        printf("Numero di SM: %d\n", device_properties.multiProcessorCount);
        printf("Numero massimo di thread per blocco: %d\n", device_properties.maxThreadsPerBlock);
        printf("Numero massimo di thread per SM: %d\n", device_properties.maxThreadsPerMultiProcessor);
        printf("Numero di CUDA cores per SM: %d\n", device_properties.maxThreadsPerMultiProcessor); // Un'accurata approssimazione
        printf("Memoria totale (in GB): %.2f GB\n", device_properties.totalGlobalMem / (1024.0 * 1024.0 * 1024.0));
    }

    return 0;
}
