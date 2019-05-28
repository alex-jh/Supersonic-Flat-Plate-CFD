#include "supersonic_flow.h"
#include "supersonic_cone.h"
#include "math.h"
#include "step_size_calculator.h"
#include <algorithm>
#include <iostream>
#include "file_writer.h"
#include <chrono>

namespace supersonic_flow {

	const double pi = acos(-1.0);

	SupersonicFlow::SupersonicFlow(int imax, int jmax) : imax_(imax), jmax_(jmax), maccormack_solver_(imax, jmax, true)
	{
		flow_parameters_.mu = 1.7894e-5;
		flow_parameters_.T_0 = 288.16;

		double M_inf = 3;
		double a_inf = 340.28;
		double u_inf = a_inf * M_inf;
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

		maxit = 40000;

		T_ = Array2D<double>(imax_, jmax_);
		u_ = Array2D<double>(imax_, jmax_);
		v_ = Array2D<double>(imax_, jmax_);
		P_ = Array2D<double>(imax_, jmax_);
		rho_ = Array2D<double>(imax_, jmax_);
		M_ = Array2D<double>(imax_, jmax_);
		e_ = Array2D<double>(imax_, jmax_);
		outside_ = Array2D<NODE_TYPE>(imax_, jmax_);

		deltax_ = CalcXStep(flow_parameters_, imax_);
		deltay_ = CalcYStep(flow_parameters_, jmax_);
	}


	SupersonicFlow::~SupersonicFlow()
	{
	}

	double SupersonicFlow::CalcXStep(const FlowParameters& params, int size) {
		return params.plate_length / size;
	}

	double SupersonicFlow::CalcYStep(const FlowParameters& params, int size) {
		double rho = params.P_inf / (params.T_inf * params.R);
		double u = params.M_inf * params.a_inf;
		double Re = rho * u * params.plate_length / params.mu;
		double delta = 5 * params.plate_length / sqrt(Re);

		double lvert = 5 * delta;
		return lvert / size;
	}

	void SupersonicFlow::Run() {

		auto start_time = std::chrono::high_resolution_clock::now();

		// Sets initial conditions.
		InitializeFlowFieldVariables();
		BoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

		int mod = 10;

		maxit = 1000;

		for (int it = 0; it < maxit; it++) {

			if (it % mod == 0)
				std::cout << "Iteration " << it << " ";

			double delta_t = CalcTStep(imax_, jmax_, deltax_, deltay_, flow_parameters_, u_, v_, rho_, P_, T_, 0.6);
			double delta_t2 = CalcTStep2(imax_, jmax_, deltax_, deltay_, flow_parameters_, u_, v_, rho_, P_, T_, 0.6);

			//std::cout << delta_t << " " << delta_t2 << std::endl;

			rho_old_ = rho_;
			u_old_ = u_;
			v_old_ = v_;
			P_old_ = P_;

			maccormack_solver_.UpdatePredictor(delta_t2, deltax_, deltay_, imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_, outside_);

			BoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

			maccormack_solver_.UpdateCorrector(delta_t2, deltax_, deltay_, imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_, outside_);

			BoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

			double diff = 0.0;
			if (CheckConvergence(diff)) {
				break;
			}

			if (it % mod == 0) {
				std::cout << diff << std::endl;

				auto end_time = std::chrono::high_resolution_clock::now();
				auto time = end_time - start_time;

				std::cout << "Iterations: " << it << std::endl;
				std::cout << "Time: " << time / std::chrono::milliseconds(1) << "ms to run.\n";

			}
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		auto time = end_time - start_time;

		std::cout << "Iterations: " << maxit << std::endl;
		std::cout << "Time: " << time / std::chrono::milliseconds(1) << "ms to run.\n";

		WriteInFile(P_, deltax_, deltay_, "Pressure");
		WriteInFile(rho_, deltax_, deltay_, "Density");
		WriteInFile(T_, deltax_, deltay_, "Temperature");
		WriteInFile(u_, deltax_, deltay_, "VelocityX");
		WriteInFile(v_, deltax_, deltay_, "VelocityY");
	}

	bool SupersonicFlow::CheckConvergence(double& diff) {

		for (int i = 0; i < imax_; i++) {
			for (int j = 0; j < jmax_; j++) {
				double tmp = abs(rho_.Get(i, j) - rho_old_.Get(i, j));
				tmp = std::max(abs(u_.Get(i, j) - u_old_.Get(i, j)), tmp);
				tmp = std::max(abs(v_.Get(i, j) - v_old_.Get(i, j)), tmp);
				tmp = std::max(abs(P_.Get(i, j) - P_old_.Get(i, j)), tmp);
				diff = std::max(diff, tmp);
			}
		}

		return diff < 1.0e-6;
	}
}

