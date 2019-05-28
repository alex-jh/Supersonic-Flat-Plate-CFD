#include "supersonic_rocket_nozzle_gpu.h"


// #include "supersonic_cone.h"
#include "math.h"
#include "step_size_calculator.h"
#include <algorithm>
#include <iostream>
#include "flow_utils.h"
// #include "file_writer.h"

namespace supersonic_rocket_nozzle_gpu {

	const int IMAX = 1000;
	const int JMAX = 250;
	const double pi = acos(-1.0);
	const double angle = 0;
	const int nozzle_exit_height = JMAX / 4;
	const int nozzle_wall_height = nozzle_exit_height / 4;
	const int freestream_inlet = JMAX - nozzle_exit_height - nozzle_wall_height;
	const int nozzle_wall_width = IMAX/4;
	const int nozzle_exit_width = nozzle_wall_width - nozzle_wall_width/4;

	SupersonicRocketNozzleGpu::SupersonicRocketNozzleGpu() : SupersonicFlowGpu(IMAX, JMAX)
	{
		T_exit_ = 3.0;
		P_exit_ = 2.5;

		double M_inf = 1.2;
		double a_inf = 340.28;
		double u_inf = a_inf * M_inf;
		double plate_length = 20.0;
		double R = 287;
		double P_inf = 101325;
		double T_inf = 288.16;
		double rho_inf = P_inf / T_inf / R;
		double mu_inf = ViscositySutherlandLaw(flow_parameters_, T_inf);
		double L = sqrt(mu_inf * plate_length / rho_inf / u_inf);
		flow_parameters_.plate_length = 10.0 / L;

		freestream_inlet_ = freestream_inlet;
		nozzle_exit_width_ = nozzle_exit_width;
		nozzle_wall_width_ = nozzle_wall_width;
		nozzle_exit_height_ = nozzle_exit_height;
		nozzle_wall_height_ = nozzle_wall_height;

		cudaMemcpy(flow_parameters_dev_, &flow_parameters_, sizeof(FlowParameters), cudaMemcpyHostToDevice);

		deltax_ = CalcXStep(flow_parameters_, imax_);
		deltay_ = CalcYStep(flow_parameters_, jmax_);

		std::cout << deltax_ << " " << deltay_ << std::endl;
		std::cout << deltax_ * imax_ << " " << deltay_ * jmax_ << std::endl;

		cudaMalloc((void **)&diff_, sizeof(double));
	}


	SupersonicRocketNozzleGpu::~SupersonicRocketNozzleGpu()
	{
		cudaFree(diff_);
		cudaFree(diff_);
	}

	double SupersonicRocketNozzleGpu::CalcXStep(const FlowParameters& params, int size) {
		std::cout << params.plate_length << " " << size << std::endl;
		double deltax = (2.5*params.plate_length) / size;

		deltax_v_->memsetToArray(0, 2* nozzle_wall_width, deltax, numThreads_, numBlocks_);
		deltax_v_->memsetToArrayLinear(2 * nozzle_wall_width, IMAX / 2, deltax, 5 * deltax, numThreads_, numBlocks_);
		deltax_v_->memsetToArrayLinear(IMAX / 2, IMAX, 5 * deltax, 10 * deltax, numThreads_, numBlocks_);

		deltax_v_->partialSum(0, IMAX, 1, 1);

		return deltax;
	}

	double SupersonicRocketNozzleGpu::CalcYStep(const FlowParameters& params, int size) {		
		double rho = params.P_inf / (params.T_inf * params.R);
		double u = params.M_inf * params.a_inf;
		double Re = rho * u * params.plate_length / params.mu;
		double delta = 5 * params.plate_length / sqrt(Re);

		double lvert = 5 * delta;
		double deltay = (2*lvert) / size;

		deltay_v_->memsetToArray(0, nozzle_exit_height*2, deltay, numThreads_, numBlocks_);
		deltay_v_->memsetToArrayLinear(nozzle_exit_height * 2, JMAX, deltay, 10 * deltay, numThreads_, numBlocks_);

		deltay_v_->partialSum(0, JMAX, 1, 1);

		return deltay;
	}

}