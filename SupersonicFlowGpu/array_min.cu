#include "array_min.cuh"
#include <math.h>
#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

template 
void reduceMin<double>(int size, int threads, int blocks, int whichKernel, const Array2DGpu<double>* dat, double *d_odata);

//template void
//reduceMin<double>(int size, int threads, int blocks, int whichKernel, int imax, int jmax, double deltax, double deltay,
//	const FlowParameters* params, const Array2DGpu<double>* u, const Array2DGpu<double>* v,
//	const Array2DGpu<double>* rho, const Array2DGpu<double>* P, const Array2DGpu<double>* T, double K, double *d_odata);

bool isPow2(unsigned int x)
{
	return ((x&(x - 1)) == 0);
}

unsigned int nextPow2(unsigned int x)
{
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory {
	__device__ inline operator T *() {
		extern __shared__ int __smem[];
		return (T *)__smem;
	}

	__device__ inline operator const T *() const {
		extern __shared__ int __smem[];
		return (T *)__smem;
	}
};

// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double> {
	__device__ inline operator double *() {
		extern __shared__ double __smem_d[];
		return (double *)__smem_d;
	}

	__device__ inline operator const double *() const {
		extern __shared__ double __smem_d[];
		return (double *)__smem_d;
	}
};

__device__ int fast_mod(const int input, const int ceil) {
	// apply the modulo operator only when needed
	// (i.e. when the input is greater than the ceiling)
	return input >= ceil ? input % ceil : input;
	// NB: the assumption here is that the numbers are positive
}

/*
 This version adds multiple elements per thread sequentially.  This reduces the overall
 cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
 (Brent's Theorem optimization)

 Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
 In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
 If blockSize > 32, allocate blockSize*sizeof(T) bytes.
 */
template<class T, unsigned int blockSize, bool nIsPow2>
__global__ void reduceMin6(T *g_odata, unsigned int n, const Array2DGpu<double>* data) {
	T *sdata = SharedMemory<T>();

	// perform first level of reduction,
	// reading from global memory, writing to shared memory
	unsigned int tid = threadIdx.x;
	unsigned int idx = blockIdx.x * blockSize * 2 + threadIdx.x;
	unsigned int gridSize = blockSize * 2 * gridDim.x;

	T myMin = 99999;

	// we reduce multiple elements per thread.  The number is determined by the
	// number of active thread blocks (via gridDim).  More blocks will result
	// in a larger gridSize and therefore fewer elements per thread
	while (idx < n) {
		myMin = MIN(data->Get(idx), myMin);

		// ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
		if (nIsPow2 || idx + blockSize < n) {
			myMin = MIN(data->Get(idx + blockSize), myMin);
		}

		idx += gridSize;
	}

	// each thread puts its local sum into shared memory
	sdata[tid] = myMin;
	__syncthreads();

	// do reduction in shared mem
	if ((blockSize >= 512) && (tid < 256)) {
		sdata[tid] = myMin = MIN(sdata[tid + 256], myMin);
	}

	__syncthreads();

	if ((blockSize >= 256) && (tid < 128)) {
		sdata[tid] = myMin = MIN(sdata[tid + 128], myMin);
	}

	__syncthreads();

	if ((blockSize >= 128) && (tid < 64)) {
		sdata[tid] = myMin = MIN(sdata[tid + 64], myMin);
	}

	__syncthreads();

#if (__CUDA_ARCH__ >= 300 )
	if (tid < 32) {
		// Fetch final intermediate sum from 2nd warp
		if (blockSize >= 64) {
			myMin = MIN(sdata[tid + 32], myMin);
		}
		// Reduce final warp using shuffle
		for (int offset = warpSize / 2; offset > 0; offset /= 2) {
			float tempMyMin = __shfl_down(myMin, offset);

			myMin = MIN(tempMyMin, myMin);
		}

	}
#else
	// fully unroll reduction within a single warp
	if ((blockSize >= 64) && (tid < 32))
	{
		sdata[tid] = myMin = MIN(sdata[tid + 32], myMin);
	}

	__syncthreads();

	if ((blockSize >= 32) && (tid < 16))
	{
		sdata[tid] = myMin = MIN(sdata[tid + 16], myMin);
	}

	__syncthreads();

	if ((blockSize >= 16) && (tid < 8))
	{
		sdata[tid] = myMin = MIN(sdata[tid + 8], myMin);
	}

	__syncthreads();

	if ((blockSize >= 8) && (tid < 4))
	{
		sdata[tid] = myMin = MIN(sdata[tid + 4], myMin);
	}

	__syncthreads();

	if ((blockSize >= 4) && (tid < 2))
	{
		sdata[tid] = myMin = MIN(sdata[tid + 2], myMin);
	}

	__syncthreads();

	if ((blockSize >= 2) && (tid < 1))
	{
		sdata[tid] = myMin = MIN(sdata[tid + 1], myMin);
	}
#endif

	__syncthreads();
	// write result for this block to global mem
	if (tid == 0) {
		g_odata[blockIdx.x] = myMin;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template<class T>
void reduceMin(int size, int threads, int blocks, int whichKernel, const Array2DGpu<double>* data, T *d_odata) {
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	int smemSize =
		(threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

	if (isPow2(size)) {
		switch (threads) {
		case 512:
			reduceMin6<T, 512, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 256:
			reduceMin6<T, 256, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 128:
			reduceMin6<T, 128, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 64:
			reduceMin6<T, 64, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 32:
			reduceMin6<T, 32, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 16:
			reduceMin6<T, 16, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 8:
			reduceMin6<T, 8, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 4:
			reduceMin6<T, 4, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 2:
			reduceMin6<T, 2, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 1:
			reduceMin6<T, 1, true> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;
		}
	}
	else {
		switch (threads) {
		case 512:
			reduceMin6<T, 512, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 256:
			reduceMin6<T, 256, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 128:
			reduceMin6<T, 128, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 64:
			reduceMin6<T, 64, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 32:
			reduceMin6<T, 32, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 16:
			reduceMin6<T, 16, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 8:
			reduceMin6<T, 8, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 4:
			reduceMin6<T, 4, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 2:
			reduceMin6<T, 2, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;

		case 1:
			reduceMin6<T, 1, false> << <dimGrid, dimBlock, smemSize >> > (d_odata,
				size, data);
			break;
		}
	}

}

ArrayMin::ArrayMin(int imax, int jmax) {
	cudaMalloc((void **)&delta_t_per_block, max_blocks * sizeof(double));
	cudaMalloc((void**)&data_, sizeof(Array2DGpu<double>));
	initArray2D<double> << <1, 1 >> > ((void*)data_, imax, jmax);
	n = imax * jmax;
}

ArrayMin::~ArrayMin() {
	cudaFree(delta_t_per_block);
	destroyArray2D<double> << <1, 1 >> > (data_);
	cudaFree(data_);
}

void ArrayMin::reduceArrayMin(int threads, int blocks, int whichKernel, double* val) {

	reduceMin(n, threads, blocks, whichKernel, data_, delta_t_per_block);

	getMin << <1, blocks / 2, 1 >> > (delta_t_per_block, blocks, val);

	//double tmp = 0;
	//cudaMemcpy(&tmp, &delta_t_per_block[0], sizeof(double), cudaMemcpyDeviceToHost);

	//printf("Delta_t: %f\n", tmp);
}

//void ArrayMin::reduceArrayMin(double* arr, int n, double* val) {
//	getMin << <1, blocks / 2, 1 >> > (arr, n, double* val);
//}

__global__ void getMin(double* arr, int n, double* val) {
	int idx = threadIdx.x;
	for (int i = n / 2; i > 0; i >>= 1) {
		if (idx < i) {
			arr[idx] = MIN(arr[idx], arr[idx + i]);
			__syncthreads();
		}
	}

	if (idx == 0) {
		*val = arr[0];
	}

	//if (threadIdx.x == 0) {
	//	for (int i = 1; i < n; i++) {
	//		arr[0] = MIN(arr[0], arr[i]);
	//	}
	//}
}
