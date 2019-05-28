#include "supersonic_flow_gpu.h"
//#include "supersonic_cone.h"
//#include "math.h"
//#include <algorithm>
#include <iostream>
#include "file_writer.h"
#include "flow_utils.h"
#include <chrono>
#include <sstream>

namespace supersonic_flow_gpu {

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

		if (whichKernel == 6) {
			blocks = MIN(maxBlocks, blocks);
		}
	}

	SupersonicFlowGpu::SupersonicFlowGpu(int imax, int jmax) : imax_(imax), jmax_(jmax), maccormack_solver_(imax, jmax, true), array_min_(imax, jmax) {

		cudaSetDevice(0);
		cudaMalloc((void**)&flow_parameters_dev_, sizeof(FlowParameters));

		flow_parameters_.mu = 1.7894e-5;
		flow_parameters_.T_0 = 288.16;

		double M_inf = 1.2;
		double a_inf = 340.28;
		double u_inf = a_inf;
		double plate_length = 0.0001;
		double R = 287;
		double P_inf = 101325;
		double T_inf = 288.16;
		double rho_inf = P_inf / T_inf / R;
		double mu_inf = ViscositySutherlandLaw(flow_parameters_, T_inf);
		double L = sqrt(mu_inf*plate_length / rho_inf / u_inf);

		flow_parameters_.M_inf = M_inf;
		flow_parameters_.plate_length = 0.00001 / L;
		flow_parameters_.a_inf = a_inf / u_inf;
		flow_parameters_.P_inf = P_inf / P_inf;
		flow_parameters_.T_inf = T_inf / T_inf;
		flow_parameters_.T_wall = flow_parameters_.T_inf;
		flow_parameters_.gamma = 1.4;
		flow_parameters_.R = 287 / u_inf / u_inf * T_inf;
		flow_parameters_.Pr = 0.71;
		flow_parameters_.cv = 0.718;
		flow_parameters_.cp = 1.01;
		flow_parameters_.mu /= (P_inf * (L / u_inf));
		flow_parameters_.T_0 /= T_inf;

		cudaMemcpy(flow_parameters_dev_, &flow_parameters_, sizeof(FlowParameters), cudaMemcpyHostToDevice);
	
		maxit_ = 40000;

		int size = imax * jmax;
		cudaMalloc((void**)&T_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&u_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&v_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&P_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&rho_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&e_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&outside_, sizeof(Array2DGpu<NODE_TYPE>));

		cudaMalloc((void**)&u_old_, sizeof(Array2DGpu<double>));
		//cudaMalloc((void**)&v_old_, sizeof(Array2DGpu<double>));
		//cudaMalloc((void**)&P_old_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&rho_old_, sizeof(Array2DGpu<double>));

		cudaMalloc((void**)&deltax_v_, sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&deltay_v_, sizeof(Array2DGpu<double>));

		Init();

		int maxThreads = 256;  // number of threads per block
		int whichKernel = 6;
		int maxBlocks = 64;

		getNumBlocksAndThreads(whichKernel, imax_*jmax_, maxBlocks, maxThreads, numBlocks_,
			numThreads_);

		maccormack_solver_.numThreads_ = numThreads_;
		maccormack_solver_.numBlocks_ = numBlocks_;

		whichKernel_ = whichKernel;
	}

	SupersonicFlowGpu::~SupersonicFlowGpu() {
		Destroy();

		cudaFree(T_);
		cudaFree(u_);
		cudaFree(v_);
		cudaFree(P_);
		cudaFree(rho_);
		cudaFree(e_);
		cudaFree(outside_);

		cudaFree(u_old_);
		//cudaFree(v_old_);
		//cudaFree(P_old_);
		cudaFree(rho_old_);

		cudaFree(deltax_v_);
		cudaFree(deltay_v_);
	}

	//void SupersonicFlowGpu::Run() {
	//	addWithCuda();
	//}

	//const double pi = acos(-1.0);

	//SupersonicFlow::SupersonicFlow(int imax, int jmax) : imax_(imax), jmax_(jmax), maccormack_solver_(imax, jmax, true)
	//{
	//	flow_parameters_.mu = 1.7894e-5;
	//	flow_parameters_.T_0 = 288.16;

	//	double M_inf = 3;
	//	double a_inf = 340.28;
	//	double u_inf = a_inf * M_inf;
	//	double plate_length = 0.0001;
	//	double R = 287;
	//	double P_inf = 101325;
	//	double T_inf = 288.16;
	//	double rho_inf = P_inf / T_inf / R;
	//	double mu_inf = ViscositySutherlandLaw(flow_parameters_, T_inf);
	//	double L = sqrt(mu_inf*plate_length / rho_inf / u_inf);

	//	flow_parameters_.M_inf = M_inf;
	//	flow_parameters_.plate_length = 0.00001 / L;
	//	flow_parameters_.a_inf = a_inf / u_inf;
	//	flow_parameters_.P_inf = P_inf / P_inf;
	//	flow_parameters_.T_inf = T_inf / T_inf;
	//	flow_parameters_.T_wall = flow_parameters_.T_inf;
	//	flow_parameters_.gamma = 1.4;
	//	flow_parameters_.R = 287 / u_inf / u_inf * T_inf;
	//	flow_parameters_.Pr = 0.71;
	//	flow_parameters_.cv = 0.718;
	//	flow_parameters_.cp = 1.01;
	//	flow_parameters_.mu /= (P_inf * (L / u_inf));
	//	flow_parameters_.T_0 /= T_inf;

	//	maxit = 40000;

	//	T_ = Array2D<double>(imax_, jmax_);
	//	u_ = Array2D<double>(imax_, jmax_);
	//	v_ = Array2D<double>(imax_, jmax_);
	//	P_ = Array2D<double>(imax_, jmax_);
	//	rho_ = Array2D<double>(imax_, jmax_);
	//	M_ = Array2D<double>(imax_, jmax_);
	//	e_ = Array2D<double>(imax_, jmax_);
	//	outside_ = Array2D<NODE_TYPE>(imax_, jmax_);

	//	deltax_ = CalcXStep(flow_parameters_, imax_);
	//	deltay_ = CalcYStep(flow_parameters_, jmax_);
	//}


	//SupersonicFlow::~SupersonicFlow()
	//{
	//}

	//double SupersonicFlowGpu::CalcXStep(const FlowParameters& params, int size) {
	//	return params.plate_length / size;
	//}

	//double SupersonicFlowGpu::CalcYStep(const FlowParameters& params, int size) {
	//	double rho = params.P_inf / (params.T_inf * params.R);
	//	double u = params.M_inf * params.a_inf;
	//	double Re = rho * u * params.plate_length / params.mu;
	//	double delta = 5 * params.plate_length / sqrt(Re);

	//	double lvert = 5 * delta;
	//	return lvert / size;
	//}

	void SupersonicFlowGpu::Run() {

		auto start_time = std::chrono::high_resolution_clock::now();

		// Sets initial conditions.
		InitializeFlowFieldVariables();
		cudaDeviceSynchronize();
		BoundaryConditions();
		cudaDeviceSynchronize();

		int mod = 10000;
		int mod_f = 100000;

		maxit_ = 1500000;

		for (int it = 0; it < maxit_; it++) {

			if (it % mod == 0)
				std::cout << "Iteration " << it << " " << std::endl;


			//double delta_t = CalcTStep(imax_, jmax_, deltax_, deltay_, flow_parameters_, u_, v_, rho_, P_, T_, 0.6);
			//double delta_t2 = CalcTStep2(imax_, jmax_, deltax_, deltay_, flow_parameters_, u_, v_, rho_, P_, T_, 0.6);
			double* delta_t2 = step_size_calculator_.CalcTStep2(imax_, jmax_, deltax_, deltay_, flow_parameters_dev_, u_, v_, rho_, P_, T_, 0.6,
				numThreads_, numBlocks_, whichKernel_, array_min_);

			//std::cout << delta_t2 << std::endl;

			//rho_old_->copyArray2D(imax_*jmax_, rho_, numThreads_, numBlocks_);
			//rho_old_ = rho_;
			u_old_->copyArray2D(imax_*jmax_, u_, numThreads_, numBlocks_);
			//v_old_ = v_;
			//P_old_ = P_;

			maccormack_solver_.UpdatePredictor(delta_t2, deltax_v_, deltay_v_, imax_, jmax_, flow_parameters_dev_, u_, v_, rho_, P_, T_, e_, outside_);

			BoundaryConditions();
			// cudaDeviceSynchronize();

			maccormack_solver_.UpdateCorrector(delta_t2, deltax_v_, deltay_v_, imax_, jmax_, flow_parameters_dev_, u_, v_, rho_, P_, T_, e_, outside_);

			BoundaryConditions();
			// cudaDeviceSynchronize();

			if (it % mod == 0) {
				auto end_time = std::chrono::high_resolution_clock::now();
				auto time = end_time - start_time;
				std::cout << "Time: " << time / std::chrono::milliseconds(1) << "ms to run.\n";
				double diff = 0.0;
				if (CheckConvergence(diff)) {
					std::cout << diff << std::endl;
					break;
				}
				std::cout << "Convergence diff: " << diff << std::endl;
			}

			if (it % mod_f == 0) {
				std::stringstream ss;
				ss << "_" << it;
				WriteToFile(ss.str());
			}

#ifdef DEBUG
			cudaDeviceSynchronize();
#endif
		}

		cudaDeviceSynchronize();

		auto end_time = std::chrono::high_resolution_clock::now();
		auto time = end_time - start_time;

		std::cout << "Iterations: " << maxit_ << std::endl;
		std::cout << "Time: " << time / std::chrono::milliseconds(1) << "ms to run.\n";

		WriteToFile("");

	}

	void SupersonicFlowGpu::WriteToFile(std::string s) {
		WriteInFile<double>(P_, imax_, jmax_, deltax_v_, deltay_v_, "Pressure" + s);
		WriteInFile<double>(rho_, imax_, jmax_, deltax_v_, deltay_v_, "Density" + s);
		WriteInFile<double>(T_, imax_, jmax_, deltax_v_, deltay_v_, "Temperature" + s);
		WriteInFile<double>(u_, imax_, jmax_, deltax_v_, deltay_v_, "VelocityX" + s);
		WriteInFile<double>(v_, imax_, jmax_, deltax_v_, deltay_v_, "VelocityY" + s);
		//WriteInFile<NODE_TYPE>(outside_, imax_, jmax_, deltax_v_, deltay_v_, "Type" + s);
	}
}
