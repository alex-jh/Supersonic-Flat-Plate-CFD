#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "array2d_gpu.cuh"
#include "flow_parameters.h"
#include "step_size_calculator.h"
#include "maccormack_solver_gpu.h"
#include <string>

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

namespace supersonic_flow_gpu {

	__global__ void addKernel(Array2DGpu<int>* c, const Array2DGpu<int>* a, const Array2DGpu<int>* b);

	__global__ void copyData(unsigned int n, int imax, int jmax, const Array2DGpu<double>* x,
		const Array2DGpu<double>* y, Array2DGpu<double>* data);

	class SupersonicFlowGpu
	{
	public:
		SupersonicFlowGpu(int imax, int jmax);
		~SupersonicFlowGpu();

		void Init();
		void Destroy();

		void Run();

		// cudaError_t addWithCuda();

		virtual void InitializeFlowFieldVariables() = 0;

		virtual void BoundaryConditions() = 0;

		bool CheckConvergence(double& diff);

		//cudaError_t SupersonicFlowGpu::addWithCuda(int *c, const int *a, const int *b, unsigned int size);

	private:
		virtual double CalcXStep(const FlowParameters& params, int size) = 0;
		virtual double CalcYStep(const FlowParameters& params, int size) = 0;

		void WriteToFile(std::string s);

	protected:
		FlowParameters flow_parameters_;
		FlowParameters* flow_parameters_dev_;

		MacCormackSolverGpu maccormack_solver_;
		int imax_;
		int jmax_;
		int maxit_;
		double deltax_;
		double deltay_;

		Array2DGpu<double>* T_;
		Array2DGpu<double>* u_;
		Array2DGpu<double>* v_;
		Array2DGpu<double>* P_;
		Array2DGpu<double>* rho_;
		Array2DGpu<double>* e_;

		//Array2D<double> T_old_;
		Array2DGpu<double>* u_old_;
		//Array2DGpu<double>* v_old_;
		//Array2DGpu<double>* P_old_;
		Array2DGpu<double>* rho_old_;
		//Array2D<double> e_old_;

		Array2DGpu<NODE_TYPE>* outside_;

		int numThreads_;
		int numBlocks_;
		int whichKernel_;

		StepSizeCalculator step_size_calculator_;
		ArrayMin array_min_;
		double* diff_;

		Array2DGpu<double>* deltax_v_;
		Array2DGpu<double>* deltay_v_;
	};
}

