#include "supersonic_flow_gpu.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace supersonic_flow_gpu {
	__global__ void addKernel(Array2DGpu<int>* c, const Array2DGpu<int>* a, const Array2DGpu<int>* b)
	{
		c->Set(threadIdx.x, blockIdx.x, a->Get(threadIdx.x, blockIdx.x) + b->Get(threadIdx.x, blockIdx.x));
	}

	void SupersonicFlowGpu::Init() {
		initArray2D<double> << <1, 1 >> > ((void*)T_, imax_, jmax_);
		initArray2D<double> << <1, 1 >> > ((void*)u_, imax_, jmax_);
		initArray2D<double> << <1, 1 >> > ((void*)v_, imax_, jmax_);
		initArray2D<double> << <1, 1 >> > ((void*)P_, imax_, jmax_);
		initArray2D<double> << <1, 1 >> > ((void*)rho_, imax_, jmax_);
		initArray2D<double> << <1, 1 >> > ((void*)e_, imax_, jmax_);
		initArray2D<NODE_TYPE> << <1, 1 >> > ((void*)outside_, imax_, jmax_);

		initArray2D<double> << <1, 1 >> > ((void*)u_old_, imax_, jmax_);
		//initArray2D<double> << <1, 1 >> > ((void*)v_old_, imax_, jmax_);
		//initArray2D<double> << <1, 1 >> > ((void*)P_old_, imax_, jmax_);
		initArray2D<double> << <1, 1 >> > ((void*)rho_old_, imax_, jmax_);

		initArray2D<double> << <1, 1 >> > ((void*)deltax_v_, 2, imax_);
		initArray2D<double> << <1, 1 >> > ((void*)deltay_v_, 2, jmax_);

		cudaDeviceSynchronize();
	}

	void SupersonicFlowGpu::Destroy() {
		destroyArray2D<double> << <1, 1 >> > (T_);
		destroyArray2D<double> << <1, 1 >> > (u_);
		destroyArray2D<double> << <1, 1 >> > (v_);
		destroyArray2D<double> << <1, 1 >> > (P_);
		destroyArray2D<double> << <1, 1 >> > (rho_);
		destroyArray2D<double> << <1, 1 >> > (e_);
		destroyArray2D<NODE_TYPE> << <1, 1 >> > (outside_);

		destroyArray2D<double> << <1, 1 >> > (u_old_);
		//destroyArray2D<double> << <1, 1 >> > (v_old_);
		//destroyArray2D<double> << <1, 1 >> > (P_old_);
		destroyArray2D<double> << <1, 1 >> > (rho_old_);

		destroyArray2D<double> << <1, 1 >> > (deltax_v_);
		destroyArray2D<double> << <1, 1 >> > (deltay_v_);

		cudaDeviceSynchronize();
	}

	//// Helper function for using CUDA to add vectors in parallel.
	//cudaError_t SupersonicFlowGpu::addWithCuda()
	//{
	//	cudaError_t cudaStatus;

	//	cudaSetDevice(0);

	//	printf("hola\n");

	//	addKernel << <imax_, jmax_ >> > (c_, a_, b_);

	//	cudaStatus = cudaGetLastError();

	//	cudaDeviceSynchronize();

	//	cudaStatus = cudaGetLastError();

	//	int* tmp;
	//	cudaMalloc((void**)&tmp, imax_ * jmax_ * sizeof(int));

	//	copyArray2D<int> << <imax_, jmax_ >> > (a_, tmp);

	//	int* tmp_local = new int[imax_*jmax_];
	//	cudaMemcpy(tmp_local, tmp, imax_*jmax_ * sizeof(int), cudaMemcpyDeviceToHost);

	//	for (int i = 0; i < imax_*jmax_; i++) printf("%d ", tmp_local[i]);

	//	int hh = 0;

	//	hh++;

	//	delete tmp;

	//	return cudaStatus;
	//}

	__global__ void copyData(unsigned int n, int imax, int jmax, const Array2DGpu<double>* x, 
		const Array2DGpu<double>* y, Array2DGpu<double>* data) {

		// perform first level of reduction,
		// reading from global memory, writing to shared memory
		unsigned int tid = threadIdx.x;
		unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned int gridSize = blockDim.x * gridDim.x;

		// we reduce multiple elements per thread.  The number is determined by the
		// number of active thread blocks (via gridDim).  More blocks will result
		// in a larger gridSize and therefore fewer elements per thread
		while (idx < n) {
			unsigned int ii = idx / jmax;
			unsigned int jj = idx - ii * jmax;
			data->Set(idx, - ABS(x->Get(idx) - y->Get(idx)));

			idx += gridSize;
		}
	}

	bool SupersonicFlowGpu::CheckConvergence(double& diff) {

		copyData << <numBlocks_, numThreads_, 1 >> > (imax_ * jmax_, imax_, jmax_, u_, u_old_, array_min_.getData());

		array_min_.reduceArrayMin(numThreads_, numBlocks_, whichKernel_, diff_);

		//reduceArrayMin(delta_t_per_block, numBlocks);

		//cudaMemcpy(out, d_out, numBlocks * sizeof(float), cudaMemcpyDeviceToHost);

		//double res = out[0];
		//for (int i = 1; i < numBlocks; i++) {
		//	res = MIN(res, out[i]);
		//}
		//cudaDeviceSynchronize();
		cudaMemcpy(&diff, diff_, sizeof(double), cudaMemcpyDeviceToHost);
		diff = -diff;

		return diff < 1.0e-6;
	}
}