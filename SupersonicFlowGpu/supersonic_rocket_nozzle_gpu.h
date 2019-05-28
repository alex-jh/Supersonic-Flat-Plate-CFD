#pragma once
#include "supersonic_flow_gpu.h"
#include "flow_parameters.h"
//#include "math_utils.h"
#include "maccormack_solver_gpu.h"
#include "array2d_gpu.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace supersonic_rocket_nozzle_gpu {

	template <unsigned int blockSize>
	__global__ void InitializeFlowFieldVariablesDevice(int n, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, int nozzle_exit_width, int nozzle_exit_height,
		int nozzle_wall_width, int nozzle_wall_height);

	template <unsigned int blockSize>
	__global__ void InletBoundaryConditions(int n, int start, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside);

	template <unsigned int blockSize>
	__global__ void NozzleExitBoundaryConditions(int n, int i, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside, double P_exit, double T_exit);

	template <unsigned int blockSize>
	__global__ void NozzleWallBoundaryConditions(int n, int start_width, int end_width, int height1, int height2, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside);

	template <unsigned int blockSize>
	__global__ void SymmetryBoundaryConditions(int n, int start, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside);

	template <unsigned int blockSize>
	__global__ void TopBoundaryConditions(int n, int imax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside);

	template <unsigned int blockSize>
	__global__ void OutletBoundaryConditions(int n, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* outside);

	void getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks,
		int maxThreads, int &blocks, int &threads);

	class SupersonicRocketNozzleGpu : public supersonic_flow_gpu::SupersonicFlowGpu
	{
	public:
		SupersonicRocketNozzleGpu();
		~SupersonicRocketNozzleGpu();

		void InitializeFlowFieldVariables();

		void BoundaryConditions();

	private:
		double CalcXStep(const FlowParameters& params, int size);
		double CalcYStep(const FlowParameters& params, int size);

	private:
		double T_exit_;
		double P_exit_;

		int freestream_inlet_;
		int nozzle_exit_width_;
		int nozzle_wall_width_;
		int nozzle_exit_height_;
		int nozzle_wall_height_;
	};
}

