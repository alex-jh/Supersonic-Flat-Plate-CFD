
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdexcept>

#include "supersonic_rocket_nozzle_gpu.h"


int main() {
	cudaError_t cudaStatus = cudaDeviceSetLimit(cudaLimitMallocHeapSize, 3 * 1024ULL * 1024ULL * 1024ULL);

	supersonic_rocket_nozzle_gpu::SupersonicRocketNozzleGpu supersonic_flow;
	supersonic_flow.Run();

 //   const int arraySize = 5;
 //   const int a[arraySize] = { 1, 2, 3, 4, 5 };
 //   const int b[arraySize] = { 10, 20, 30, 40, 50 };
 //   int c[arraySize] = { 0 };

	//int maxThreads = 256;  // number of threads per block
	//int whichKernel = 6;
	//int maxBlocks = 64;

	//int testIterations = 100;
	//int num_els = 1024;


	//float* d_in = NULL;
	//float* d_out = NULL;
	//int *d_idxs = NULL;
	//int *d_oIdxs = NULL;

	//printf("%d elements\n", num_els);
	//printf("%d threads (max)\n", maxThreads);

	//int numBlocks = 0;
	//int numThreads = 0;
	//getNumBlocksAndThreads(whichKernel, num_els, maxBlocks, maxThreads, numBlocks,
	//	numThreads);

 //   // Add vectors in parallel.
 //   cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
 //   if (cudaStatus != cudaSuccess) {
 //       fprintf(stderr, "addWithCuda failed!");
 //       return 1;
 //   }

 //   printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
 //       c[0], c[1], c[2], c[3], c[4]);

 //   // cudaDeviceReset must be called before exiting in order for profiling and
 //   // tracing tools such as Nsight and Visual Profiler to show complete traces.
 //   cudaStatus = cudaDeviceReset();
 //   if (cudaStatus != cudaSuccess) {
 //       fprintf(stderr, "cudaDeviceReset failed!");
 //       return 1;
 //   }

    return 0;
}
