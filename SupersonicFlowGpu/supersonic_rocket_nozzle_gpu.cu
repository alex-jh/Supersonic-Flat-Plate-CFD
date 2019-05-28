#include "supersonic_rocket_nozzle_gpu.h"
#include "math.h"

namespace supersonic_rocket_nozzle_gpu {

	//#define PRINT_CUDA_PERF

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

	////////////////////////////////////////////////////////////////////////////////
	// Compute the number of threads and blocks to use for the given reduction kernel
	// For the kernels >= 3, we set threads / block to the minimum of maxThreads and
	// n/2. For kernels < 3, we set to the minimum of maxThreads and n.  For kernel
	// 6, we observe the maximum specified number of blocks, because each thread in
	// that kernel can process a variable number of elements.
	////////////////////////////////////////////////////////////////////////////////
	void getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks,
		int maxThreads, int &blocks, int &threads) {

		//get device capability, to avoid block/grid size exceed the upper bound
		cudaDeviceProp prop;
		int device;
		cudaGetDevice(&device);
		cudaGetDeviceProperties(&prop, device);

		if (whichKernel < 3) {
			threads = (n < maxThreads) ? nextPow2(n) : maxThreads;
			blocks = (n + threads - 1) / threads;
		}
		else {
			threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
			blocks = (n + (threads * 2 - 1)) / (threads * 2);
		}

		if ((float)threads * blocks
			> (float) prop.maxGridSize[0] * prop.maxThreadsPerBlock) {
			printf("n is too large, please choose a smaller number!\n");
		}

		if (blocks > prop.maxGridSize[0]) {
			printf(
				"Grid size <%d> exceeds the device capability <%d>, set block size as %d (original %d)\n",
				blocks, prop.maxGridSize[0], threads * 2, threads);

			blocks /= 2;
			threads *= 2;
		}

		if (whichKernel == 6 && maxBlocks < blocks) {
			blocks = maxBlocks;
		}
	}


	template <unsigned int blockSize>
	__global__ void InitializeFlowFieldVariablesDevice(int n, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, int nozzle_exit_width, int nozzle_exit_height,
		int nozzle_wall_width, int nozzle_wall_height)
	{
		unsigned int idx = blockIdx.x * blockSize + threadIdx.x;
		unsigned int gridSize = blockSize * gridDim.x;

		while (idx < n) {
			unsigned int i = idx / jmax;
			unsigned int j = idx - i * jmax;

			T->Get(i, j) = params->T_inf;
			P->Get(i, j) = params->P_inf;
			rho->Get(i, j) = params->P_inf / params->T_inf / params->R;
			u->Get(i, j) = 0;
			v->Get(i, j) = 0;
			e->Get(i, j) = params->cv * params->T_inf;

			type->Set(i, j, NODE_TYPE::INSIDE);

			if (i < nozzle_exit_width && j <= nozzle_exit_height) {
				type->Set(i, j, NODE_TYPE::OUTSIDE);
			}

			if (j < nozzle_exit_height + nozzle_wall_height && j > nozzle_exit_height &&
				i < nozzle_wall_width) {
				type->Set(i, j, NODE_TYPE::OUTSIDE);
			}

			idx += gridSize;
		}
	}

	void SupersonicRocketNozzleGpu::InitializeFlowFieldVariables() {

		dim3 dimBlock(numThreads_, 1, 1);
		dim3 dimGrid(numBlocks_, 1, 1);

#ifdef PRINT_CUDA_PERF
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start);
#endif

		switch (numThreads_) {
		case 512:
			InitializeFlowFieldVariablesDevice<512> << <dimGrid, dimBlock, 1 >> > (imax_*jmax_, jmax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_, nozzle_exit_width_, nozzle_exit_height_,
				nozzle_wall_width_, nozzle_wall_height_);
			break;

		case 256:
			InitializeFlowFieldVariablesDevice<256> << <dimGrid, dimBlock, 1 >> > (imax_*jmax_, jmax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_, nozzle_exit_width_, nozzle_exit_height_,
				nozzle_wall_width_, nozzle_wall_height_);
			break;

		case 128:
			InitializeFlowFieldVariablesDevice<128> << <dimGrid, dimBlock, 1 >> > (imax_*jmax_, jmax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_, nozzle_exit_width_, nozzle_exit_height_,
				nozzle_wall_width_, nozzle_wall_height_);
			break;

		default:
			printf("Unsupported number of threads\n\n");
		}

		//cudaDeviceSynchronize();
#ifdef PRINT_CUDA_PERF
		cudaEventRecord(stop);

		cudaEventSynchronize(stop);
		float milliseconds = 0;
		cudaEventElapsedTime(&milliseconds, start, stop);

		printf("Time: %f\n\n", milliseconds);
		printf("Effective Bandwidth (GB/s): %f\n\n", imax_*jmax_ * 4 * 3 / milliseconds / 1e6);
#endif

	}

	void SupersonicRocketNozzleGpu::BoundaryConditions() {

		dim3 dimBlock(numThreads_, 1, 1);
		dim3 dimGrid(numBlocks_, 1, 1);

		switch (numThreads_) {
		case 512:
			InletBoundaryConditions<512> << <dimGrid, dimBlock, 1 >> > (jmax_, jmax_ - freestream_inlet_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			NozzleExitBoundaryConditions<512> << <dimGrid, dimBlock, 1 >> > (jmax_, nozzle_exit_width_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_, P_exit_, T_exit_);

			NozzleWallBoundaryConditions<512> << <dimGrid, dimBlock, 1 >> > (imax_, nozzle_exit_width_, nozzle_wall_width_, 
				nozzle_exit_height_, nozzle_wall_height_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			SymmetryBoundaryConditions<512> << <dimGrid, dimBlock, 1 >> > (imax_, nozzle_exit_width_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			OutletBoundaryConditions<512> << <dimGrid, dimBlock, 1 >> > (jmax_, imax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			TopBoundaryConditions<512> << <dimGrid, dimBlock, 1 >> > (imax_, jmax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			break;

		case 256:
			InletBoundaryConditions<256> << <dimGrid, dimBlock, 1 >> > (jmax_, jmax_ - freestream_inlet_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			NozzleExitBoundaryConditions<256> << <dimGrid, dimBlock, 1 >> > (nozzle_exit_height_, nozzle_exit_width_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_, P_exit_, T_exit_);

			NozzleWallBoundaryConditions<256> << <dimGrid, dimBlock, 1 >> > (imax_, nozzle_exit_width_, nozzle_wall_width_,
				nozzle_exit_height_, nozzle_wall_height_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			SymmetryBoundaryConditions<256> << <dimGrid, dimBlock, 1 >> > (imax_, nozzle_exit_width_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			OutletBoundaryConditions<256> << <dimGrid, dimBlock, 1 >> > (jmax_, imax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			TopBoundaryConditions<256> << <dimGrid, dimBlock, 1 >> > (imax_, jmax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			break;

		case 128:
			InletBoundaryConditions<128> << <dimGrid, dimBlock, 1 >> > (jmax_, jmax_ - freestream_inlet_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			NozzleExitBoundaryConditions<128> << <dimGrid, dimBlock, 1 >> > (jmax_, nozzle_exit_width_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_, P_exit_, T_exit_);

			NozzleWallBoundaryConditions<128> << <dimGrid, dimBlock, 1 >> > (imax_, nozzle_exit_width_, nozzle_wall_width_,
				nozzle_exit_height_, nozzle_wall_height_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			SymmetryBoundaryConditions<128> << <dimGrid, dimBlock, 1 >> > (imax_, nozzle_exit_width_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			OutletBoundaryConditions<128> << <dimGrid, dimBlock, 1 >> > (jmax_, imax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			TopBoundaryConditions<128> << <dimGrid, dimBlock, 1 >> > (imax_, jmax_, flow_parameters_dev_,
				u_, v_, rho_, P_, T_, e_, outside_);

			break;

		default:
			printf("Unsupported number of threads %d.\n\n", numThreads_);
		}

	}

	template <unsigned int blockSize>
	__global__ void InletBoundaryConditions(int n, int start, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside) {

		unsigned int idx = blockIdx.x * blockSize + threadIdx.x;
		unsigned int gridSize = blockSize * gridDim.x;

		while (idx < n) {

			if (idx >= start) {
				u->Set(0, idx, 0);
				v->Set(0, idx, 0);
				P->Set(0, idx, params->P_inf);
				T->Set(0, idx, params->T_inf);
				rho->Set(0, idx, params->P_inf / params->R / params->T_inf);
				e->Set(0, idx, params->T_inf * params->cv);

				outside->Set(0, idx, NODE_TYPE::BOUNDARY);
			}

			idx += gridSize;
		}
	}

	template <unsigned int blockSize>
	__global__ void NozzleExitBoundaryConditions(int n, int i, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside, double P_exit, double T_exit) {

		unsigned int idx = blockIdx.x * blockSize + threadIdx.x;
		unsigned int gridSize = blockSize * gridDim.x;

		while (idx < n) {
			u->Set(i, idx, params->a_inf * params->M_inf);
			v->Set(i, idx, 0);
			P->Set(i, idx, P_exit);
			T->Set(i, idx, T_exit);
			rho->Set(i, idx, P_exit / params->R / T_exit);
			e->Set(i, idx, T_exit * params->cv);

			outside->Set(i, idx, NODE_TYPE::BOUNDARY);

			idx += gridSize;
		}
	}

	template <unsigned int blockSize>
	__global__ void NozzleWallBoundaryConditions(int n, int start_width, int end_width, int height1, int height2, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside) {

		unsigned int idx = blockIdx.x * blockSize + threadIdx.x;
		unsigned int gridSize = blockSize * gridDim.x;

		while (idx < end_width) {
			if (idx >= start_width) {
				double p_ij = 2 * P->Get(idx, height1 - 1) - P->Get(idx, height1 - 2);

				u->Set(idx, height1, 0);
				v->Set(idx, height1, 0);
				P->Set(idx, height1, p_ij);
				T->Set(idx, height1, params->T_wall);
				rho->Set(idx, height1, p_ij / params->T_wall / params->R);
				e->Set(idx, height1, params->cv * params->T_wall);

				outside->Set(idx, height1, NODE_TYPE::BOUNDARY);
			}

			double p_ij = 2 * P->Get(idx, height1 + height2 + 1) - P->Get(idx, height1 + height2 + 2);

			u->Set(idx, height1 + height2, 0);
			v->Set(idx, height1 + height2, 0);
			P->Set(idx, height1 + height2, p_ij);
			T->Set(idx, height1 + height2, params->T_wall);
			rho->Set(idx, height1 + height2, p_ij / params->T_wall / params->R);
			e->Set(idx, height1 + height2, params->cv * params->T_wall);

			outside->Set(idx, height1 + height2, NODE_TYPE::BOUNDARY);

			idx += gridSize;
		}

		idx = threadIdx.x * gridDim.x + blockIdx.x;

		while (idx < height1 + height2) {
			if (idx >= height1) {
				double p_ij = 2 * P->Get(end_width + 1, idx) - P->Get(end_width + 2, idx);

				u->Set(end_width, idx, 0);
				v->Set(end_width, idx, 0);
				P->Set(end_width, idx, p_ij);
				T->Set(end_width, idx, params->T_wall);
				rho->Set(end_width, idx, p_ij / params->T_wall / params->R);
				e->Set(end_width, idx, params->cv * params->T_wall);

				outside->Set(end_width, idx, NODE_TYPE::BOUNDARY);
			}

			idx += gridSize;
		}
	}

	template <unsigned int blockSize>
	__global__ void SymmetryBoundaryConditions(int n, int start, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside) {

		unsigned int idx = blockIdx.x * blockSize + threadIdx.x;
		unsigned int gridSize = blockSize * gridDim.x;

		while (idx < n) {
			if (idx > start) {
				double K = 2;

				u->Set(idx, 0, (K * K * u->Get(idx, 1) - u->Get(idx, 2)) / (K * K - 1));
				v->Set(idx, 0, 0);
				rho->Set(idx, 0, (K * K * rho->Get(idx, 1) - rho->Get(idx, 2)) / (K * K - 1));

				// c_p*T + u*u/2 is constant
				T->Set(idx, 0, (params->cp*T->Get(start, 0) + u->Get(start, 0)* u->Get(start, 0) / 2 -
					u->Get(idx, 0)* u->Get(idx, 0) / 2) / params->cp);

				//T->Set(idx, 0, (K * K * T->Get(idx, 1) - T->Get(idx, 2)) / (K * K - 1));
				P->Set(idx, 0, T->Get(idx, 0) * params->R * rho->Get(idx, 0));

				//P->Set(idx, 0, 2 * P->Get(idx, 1) - P->Get(idx, 2));
				//T->Set(idx, 0, P->Get(idx, 0) / params->R / rho->Get(idx, 0));

				e->Set(idx, 0, T->Get(idx, 0) * params->cv);

				outside->Set(idx, 0, NODE_TYPE::BOUNDARY);
			}

			idx += gridSize;
		}
	}

	template <unsigned int blockSize>
	__global__ void OutletBoundaryConditions(int n, int imax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside) {

		unsigned int idx = blockIdx.x * blockSize + threadIdx.x;
		unsigned int gridSize = blockSize * gridDim.x;

		while (idx < n) {

			double c2 = params->gamma * P->Get(imax - 2, idx) / rho->Get(imax - 2, idx);
			double vel2 = u->Get(imax - 2, idx)*u->Get(imax - 2, idx) + v->Get(imax - 2, idx)*v->Get(imax - 2, idx);

			if (vel2 > c2) {
				u->Set(imax - 1, idx, 3 * u->Get(imax - 2, idx) - 3 * u->Get(imax - 3, idx) + u->Get(imax - 4, idx));
				v->Set(imax - 1, idx, 3 * v->Get(imax - 2, idx) - 3 * v->Get(imax - 3, idx) + v->Get(imax - 4, idx));
				T->Set(imax - 1, idx, 3 * T->Get(imax - 2, idx) - 3 * T->Get(imax - 3, idx) + T->Get(imax - 4, idx));
				rho->Set(imax - 1, idx, 3 * rho->Get(imax - 2, idx) - 3 * rho->Get(imax - 3, idx) + rho->Get(imax - 4, idx));
				P->Set(imax - 1, idx, 3 * P->Get(imax - 2, idx) - 3 * P->Get(imax - 3, idx) + P->Get(imax - 4, idx));
				e->Set(imax - 1, idx, 3 * e->Get(imax - 2, idx) - 3 * e->Get(imax - 3, idx) + e->Get(imax - 4, idx));
			}
			else {
				u->Set(imax - 1, idx, u->Get(imax - 2, idx));
				v->Set(imax - 1, idx, v->Get(imax - 2, idx));
				T->Set(imax - 1, idx, T->Get(imax - 2, idx));
				rho->Set(imax - 1, idx, rho->Get(imax - 2, idx));
				P->Set(imax - 1, idx, P->Get(imax - 2, idx));
				e->Set(imax - 1, idx, e->Get(imax - 2, idx));
			}

			//u->Set(imax - 1, idx, 2 * u->Get(imax - 2, idx) - u->Get(imax - 3, idx));
			//v->Set(imax - 1, idx, 2 * v->Get(imax - 2, idx) - v->Get(imax - 3, idx));
			//T->Set(imax - 1, idx, 2 * T->Get(imax - 2, idx) - T->Get(imax - 3, idx));
			//rho->Set(imax - 1, idx, 2 * rho->Get(imax - 2, idx) - rho->Get(imax - 3, idx));
			//P->Set(imax - 1, idx, 2 * P->Get(imax - 2, idx) - P->Get(imax - 3, idx));
			//e->Set(imax - 1, idx, 2 * e->Get(imax - 2, idx) - e->Get(imax - 3, idx));			//u->Set(imax - 1, idx, 2 * u->Get(imax - 2, idx) - u->Get(imax - 3, idx));
			//v->Set(imax - 1, idx, 2 * v->Get(imax - 2, idx) - v->Get(imax - 3, idx));
			//T->Set(imax - 1, idx, 2 * T->Get(imax - 2, idx) - T->Get(imax - 3, idx));
			//rho->Set(imax - 1, idx, 2 * rho->Get(imax - 2, idx) - rho->Get(imax - 3, idx));
			//P->Set(imax - 1, idx, 2 * P->Get(imax - 2, idx) - P->Get(imax - 3, idx));
			//e->Set(imax - 1, idx, 2 * e->Get(imax - 2, idx) - e->Get(imax - 3, idx));

			outside->Set(imax - 1, idx, NODE_TYPE::BOUNDARY);

			idx += gridSize;
		}
	}

	template <unsigned int blockSize>
	__global__ void TopBoundaryConditions(int n, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside) {

		unsigned int idx = blockIdx.x * blockSize + threadIdx.x;
		unsigned int gridSize = blockSize * gridDim.x;

		while (idx < n) {

			double c2 = params->gamma * P->Get(idx, jmax - 2) / rho->Get(idx, jmax - 2);
			double vel2 = u->Get(idx, jmax - 2)*u->Get(idx, jmax - 2) + v->Get(idx, jmax - 2)*v->Get(idx, jmax - 2);

			if (vel2 > c2) {
				u->Set(idx, jmax - 1, 3 * u->Get(idx, jmax - 2) - 3 * u->Get(idx, jmax - 3) + u->Get(idx, jmax - 4));
				v->Set(idx, jmax - 1, 3 * v->Get(idx, jmax - 2) - 3 * v->Get(idx, jmax - 3) + v->Get(idx, jmax - 4));
				T->Set(idx, jmax - 1, 3 * T->Get(idx, jmax - 2) - 3 * T->Get(idx, jmax - 3) + T->Get(idx, jmax - 4));
				rho->Set(idx, jmax - 1, 3 * rho->Get(idx, jmax - 2) - 3 * rho->Get(idx, jmax - 3) + rho->Get(idx, jmax - 4));
				P->Set(idx, jmax - 1, 3 * P->Get(idx, jmax - 2) - 3 * P->Get(idx, jmax - 3) + P->Get(idx, jmax - 4));
				e->Set(idx, jmax - 1, 3 * e->Get(idx, jmax - 2) - 3 * e->Get(idx, jmax - 3) + e->Get(idx, jmax - 4));
			} else {
				u->Set(idx, jmax - 1, u->Get(idx, jmax - 2));
				v->Set(idx, jmax - 1, v->Get(idx, jmax - 2));
				T->Set(idx, jmax - 1, T->Get(idx, jmax - 2));
				rho->Set(idx, jmax - 1, rho->Get(idx, jmax - 2));
				P->Set(idx, jmax - 1, P->Get(idx, jmax - 2));
				e->Set(idx, jmax - 1, e->Get(idx, jmax - 2));
			}

			//u->Set(idx, jmax - 1, 2 * u->Get(idx, jmax - 2) - u->Get(idx, jmax - 3));
			//v->Set(idx, jmax - 1, 2 * v->Get(idx, jmax - 2) - v->Get(idx, jmax - 3));
			//T->Set(idx, jmax - 1, 2 * T->Get(idx, jmax - 2) - T->Get(idx, jmax - 3));
			//rho->Set(idx, jmax - 1, 2 * rho->Get(idx, jmax - 2) - rho->Get(idx, jmax - 3));
			//P->Set(idx, jmax - 1, 2 * P->Get(idx, jmax - 2) - P->Get(idx, jmax - 3));
			//e->Set(idx, jmax - 1, 2 * e->Get(idx, jmax - 2) - e->Get(idx, jmax - 3));

			outside->Set(idx, jmax - 1, NODE_TYPE::BOUNDARY);

			idx += gridSize;
		}

	}
}