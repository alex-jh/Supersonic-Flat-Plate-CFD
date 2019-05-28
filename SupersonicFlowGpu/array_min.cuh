#pragma once
#include "array2d_gpu.cuh"

class FlowParameters;

#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x,y) ((x > y) ? x : y)
#endif

#if (__CUDA_ARCH__ < 200)
#define int_mult(x,y)   __mul24(x,y)
#else
#define int_mult(x,y)   x*y
#endif

#define inf 0x7f800000

bool isPow2(unsigned int x);

unsigned int nextPow2(unsigned int x);

template<class T, unsigned int blockSize, bool nIsPow2>
__global__ void reduceMin6(T *g_odata, unsigned int n, const Array2DGpu<double>* data);

__device__ int fast_mod(const int input, const int ceil);

template<class T>
void reduceMin(int size, int threads, int blocks, int whichKernel, const Array2DGpu<double>* dat, T *d_odata);

__global__ void getMin(double* arr, int n, double* val);

class ArrayMin {
public:
	ArrayMin(int imax, int jmax);
	~ArrayMin();

	// void reduceArrayMin(double* arr, int n, double* val);
	void reduceArrayMin(int threads, int blocks, int whichKernel, double* val);

	double* delta_t_per_block;
	double* delta_t_per_block_local;

	const int max_blocks = 512;

	Array2DGpu<double>* getData() { return data_; }

private:
	Array2DGpu<double>* data_;
	int n;
};