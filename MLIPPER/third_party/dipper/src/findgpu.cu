#include <iostream>

int main (){
    // A simple code to detect cuda-capable GPUs
    int device;
    cudaDeviceProp prop;
    if (cudaGetDevice(&device) != cudaSuccess) return 1;
    if (cudaGetDeviceProperties(&prop, device) != cudaSuccess) return 1;
    std::cout << prop.major << prop.minor;
    return 0;
}