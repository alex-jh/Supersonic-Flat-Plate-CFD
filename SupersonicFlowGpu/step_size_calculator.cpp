#include "step_size_calculator.h"
#include "flow_parameters.h"
#include <math.h>
#include "flow_utils.h"
#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <algorithm>


////////////////////////////////////////////////////////////////////////////////
//! Compute sum reduction on CPU
//! We use Kahan summation for an accurate sum of large arrays.
//! http://en.wikipedia.org/wiki/Kahan_summation_algorithm
//!
//! @param data       pointer to input data
//! @param size       number of input data elements
////////////////////////////////////////////////////////////////////////////////
template<class T>
void reduceMINCPU(T *data, int size, T *min)
{
	*min = data[0];
	T c = (T)0.0;

	for (int i = 1; i < size; i++)
	{
		T y = data[i];
		T t = MIN(*min, y);
		(*min) = t;
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//! Compute sum reduction on CPU
//! We use Kahan summation for an accurate sum of large arrays.
//! http://en.wikipedia.org/wiki/Kahan_summation_algorithm
//!
//! @param data       pointer to input data
//! @param size       number of input data elements
////////////////////////////////////////////////////////////////////////////////
template<class T>
T reduceCPU(T *data, int size)
{
	T sum = data[0];
	T c = (T)0.0;

	for (int i = 1; i < size; i++)
	{
		T y = data[i] - c;
		T t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	return sum;
}

//// TODO: Measure speed.
//unsigned long long int my_min_max_test(int num_els) {
//
//	// timers
//
//	unsigned long long int start;
//	unsigned long long int delta = 0.0;
//
//	int maxThreads = 256;  // number of threads per block
//	int whichKernel = 6;
//	int maxBlocks = 64;
//
//	int testIterations = 100;
//
//
//
//	double* d_in = NULL;
//	double* d_out = NULL;
//
//	printf("%d elements\n", num_els);
//	printf("%d threads (max)\n", maxThreads);
//
//	int numBlocks = 0;
//	int numThreads = 0;
//	getNumBlocksAndThreads(whichKernel, num_els, maxBlocks, maxThreads, numBlocks,
//		numThreads);
//
//
//
//	//  in[1024] = 34.0f;
//	//  in[333] = 55.0f;
//	//  in[23523] = -42.0f;
//
////  cudaMalloc((void**) &d_in, size);
////  cudaMalloc((void**) &d_out, size);
////  cudaMalloc((void**) &d_idxs, num_els * sizeof(int));
//
//	cudaMalloc((void **)&d_in, num_els * sizeof(double));
//	cudaMalloc((void **)&d_out, numBlocks * sizeof(double));
//
//	double* in = (double*)malloc(num_els * sizeof(double));
//	double* out = (double*)malloc(numBlocks * sizeof(double));
//
//	for (int i = 0; i < num_els; i++) {
//		in[i] = ((double)rand() + rand() + 1) / (double)RAND_MAX;
//	}
//
//	//for (int i = 0; i < num_els; i++) {
//
//	//	printf("\n %.4f", in[i]);
//	//}
//
//	// copy data directly to device memory
//	cudaMemcpy(d_in, in, num_els * sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_out, out, numBlocks * sizeof(double), cudaMemcpyHostToDevice);
//
//	// warm-up
////  reduce<float>(num_els, numThreads, numBlocks, whichKernel, d_in, d_out);
////
////  cudaMemcpy(out, d_out, numBlocks * sizeof(float), cudaMemcpyDeviceToHost);
////
////  for(int i=0; i< numBlocks; i++)
////  printf("\nFinal Result[BLK:%d]: %f", i, out[i]);
//
////  printf("\n Reduce CPU : %f", reduceCPU<float>(in, num_els));
//
//	reduceMin<double>(num_els, numThreads, numBlocks, whichKernel, d_in, d_out);
//
//	cudaMemcpy(out, d_out, numBlocks * sizeof(double), cudaMemcpyDeviceToHost);
//
//	for (int i = 0; i < numBlocks; i++) {
//		printf("\n Reduce MIN GPU value: %f", out[i]);
//	}
//
//	double minGPU;
//	reduceMINCPU<double>(out, numBlocks, &minGPU);
//
//	printf("\n\n Reduce MIN GPU value: %f", minGPU);
//
//	double min;
//	reduceMINCPU<double>(in, num_els, &min);
//
//	printf("\n\n Reduce MIN CPU value: %f", min);
//
//	cudaFree(d_in);
//	cudaFree(d_out);
//
//	free(in);
//	free(out);
//
//	//system("pause");
//
//	return delta;
//
//}

StepSizeCalculator::StepSizeCalculator() {
	cudaMalloc((void **)&delta_t_per_block, max_blocks * sizeof(double));
	cudaMalloc((void **)&delta_t, sizeof(double));
	//double* delta_t_per_block_local = (double*)malloc(max_blocks * sizeof(double));
}

StepSizeCalculator::~StepSizeCalculator() {
	cudaFree(delta_t_per_block);
	cudaFree(delta_t);
}

//double* StepSizeCalculator::CalcTStep2(const ArrayMin& array_min, int imax, int jmax, double deltax, double deltay,
//	const FlowParameters* params, const Array2DGpu<double>* u, const Array2DGpu<double>* v,
//	const Array2DGpu<double>* rho, const Array2DGpu<double>* P, const Array2DGpu<double>* T,
//	double K, int numThreads, int numBlocks, int whichKernel) {
//
//	//std::fill(out, out + numBlocks, 1.0);
//
//	//cudaMemcpy(delta_t_per_block, out, numBlocks * sizeof(double), cudaMemcpyHostToDevice);
//
//	copyData(imax*jmax, numThreads, numBlocks, whichKernel, imax, jmax, deltax, deltay,
//		params, u, v, rho, P, T, K, array_min->getData());
//
//	//reduceMin<double>(imax*jmax, numThreads, numBlocks, whichKernel, imax, jmax, deltax, deltay,
//	//	params, u, v, rho, P, T, K, delta_t_per_block);
//
//	//reduceArrayMin(delta_t_per_block, numBlocks);
//
//	//cudaMemcpy(out, d_out, numBlocks * sizeof(float), cudaMemcpyDeviceToHost);
//
//	//double res = out[0];
//	//for (int i = 1; i < numBlocks; i++) {
//	//	res = MIN(res, out[i]);
//	//}
//
//#ifdef DEBUG
//	double tmp = 0;
//	cudaMemcpy(&tmp, &delta_t_per_block[0], sizeof(double), cudaMemcpyDeviceToHost);
//
//	printf("Delta_t: %f\n", tmp);
//#endif
//
//	return delta_t_per_block;
//}

//double CalcTStep2(int imax, int jmax, double deltax, double deltay,
//	const FlowParameters& params, const Array2D<double>& u, const Array2D<double>& v,
//	const Array2D<double>& rho, const Array2D<double>& P, const Array2D<double>& T, double K) {
//
//	double delta_t = 1.0;
//
//	for (int i = 0; i < imax; i++) {
//		for (int j = 0; j < jmax; j++) {
//
//			double u_ij = u.Get(i, j);
//			double v_ij = v.Get(i, j);
//			double a_ij = sqrt(std::max(0.0, P.Get(i, j)*params.gamma / rho.Get(i, j)));
//			double mu_ij = ViscositySutherlandLaw(params, T.Get(i, j));
//			double vv_ij = abs(4. / 3.*mu_ij*(params.gamma*mu_ij / params.Pr) / rho.Get(i, j));
//
//			double delta_t_ij = 0.0;
//
//			double x = u_ij + a_ij + 1.0 / rho.Get(i, j)*(2 * params.gamma / deltax *
//				(params.gamma*mu_ij / params.Pr) + 1.0 / deltay * sqrt(2. / 3 * mu_ij * mu_ij));
//
//			double y = v_ij + a_ij + 1.0 / rho.Get(i, j)*(2 * params.gamma / deltay *
//				(params.gamma*mu_ij / params.Pr) + 1.0 / deltax * sqrt(2. / 3 * mu_ij * mu_ij));
//
//			delta_t_ij = std::min(deltax / x, deltay / y);
//			delta_t = std::min(delta_t, delta_t_ij*K);
//		}
//	}
//
//	return delta_t;
//}
