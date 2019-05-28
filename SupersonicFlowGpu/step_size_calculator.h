#pragma once
#include "array2d_gpu.cuh"
#include "array_min.cuh"

//#define DEBUG

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

//#define inf 0x7f800000
//
//bool isPow2(unsigned int x);
//
//unsigned int nextPow2(unsigned int x);
//
//template<class T, unsigned int blockSize, bool nIsPow2>
//__global__ void reduceMin6(T *g_odata, unsigned int n, int imax, int jmax, double deltax, double deltay,
//	const FlowParameters* params, const Array2DGpu<double>* u, const Array2DGpu<double>* v,
//	const Array2DGpu<double>* rho, const Array2DGpu<double>* P, const Array2DGpu<double>* t, double K);
//
//__device__ int fast_mod(const int input, const int ceil);
//
//template<class T>
//void reduceMin(int size, int threads, int blocks, int whichKernel, int imax, int jmax, double deltax, double deltay,
//	const FlowParameters* params, const Array2DGpu<double>* u, const Array2DGpu<double>* v,
//	const Array2DGpu<double>* rho, const Array2DGpu<double>* P, const Array2DGpu<double>* t, double K, T *d_odata);

__device__ double CalcLocalTStep2(int i, int j, double deltax, double deltay,
	const FlowParameters* params, const Array2DGpu<double>* u, const Array2DGpu<double>* v,
	const Array2DGpu<double>* rho, const Array2DGpu<double>* P, const Array2DGpu<double>* T, double K);

//__global__ void getMin(double* arr, int n);

__global__ void copyData(unsigned int n, int imax, int jmax, double deltax, double deltay,
	const FlowParameters* params, const Array2DGpu<double>* u, const Array2DGpu<double>* v,
	const Array2DGpu<double>* rho, const Array2DGpu<double>* P, const Array2DGpu<double>* t, double K, Array2DGpu<double>* data);

class StepSizeCalculator {
public:
	StepSizeCalculator();
	~StepSizeCalculator();

	double* CalcTStep2(int imax, int jmax, double deltax, double deltay,
		const FlowParameters* params, const Array2DGpu<double>* u, const Array2DGpu<double>* v,
		const Array2DGpu<double>* rho, const Array2DGpu<double>* P, const Array2DGpu<double>* T,
		double K, int numThreads, int numBlocks, int whichKernel, ArrayMin& array_min);

	void reduceArrayMin(double* arr, int n);

	double* delta_t_per_block;
	double* delta_t_per_block_local;
	double* delta_t;

	const int max_blocks = 512;
};

//unsigned long long int my_min_max_test(int num_els);