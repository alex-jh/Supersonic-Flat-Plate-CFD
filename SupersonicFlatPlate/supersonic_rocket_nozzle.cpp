#include "supersonic_rocket_nozzle.h"


#include "supersonic_cone.h"
#include "math.h"
#include "step_size_calculator.h"
#include <algorithm>
#include <iostream>
#include "file_writer.h"

namespace supersonic_rocket_nozzle {

	const int IMAX = 200;
	const int JMAX = 80;
	const double pi = acos(-1.0);
	const double angle = 0;
	const int nozzle_exit_height = JMAX / 4;
	const int nozzle_wall_height = nozzle_exit_height / 4;
	const int freestream_inlet = JMAX - nozzle_exit_height - nozzle_wall_height;
	const int nozzle_wall_width = IMAX / 40 + 20;
	const int nozzle_exit_width = nozzle_wall_width - 5;

	SupersonicRocketNozzle::SupersonicRocketNozzle() : SupersonicFlow(IMAX, JMAX)
	{
		T_exit_ = 2.0;
		P_exit_ = 1.5;

		double M_inf = 3;
		double a_inf = 340.28;
		double u_inf = a_inf * M_inf;
		double plate_length = 5.0;
		double R = 287;
		double P_inf = 101325;
		double T_inf = 288.16;
		double rho_inf = P_inf / T_inf / R;
		double mu_inf = ViscositySutherlandLaw(flow_parameters_, T_inf);
		double L = sqrt(mu_inf * plate_length / rho_inf / u_inf);
		flow_parameters_.plate_length = 10.0 / L;
	}


	SupersonicRocketNozzle::~SupersonicRocketNozzle()
	{
	}

	double SupersonicRocketNozzle::CalcXStep(const FlowParameters& params, int size) {
		return params.plate_length / size;
	}

	double SupersonicRocketNozzle::CalcYStep(const FlowParameters& params, int size) {
		double rho = params.P_inf / (params.T_inf * params.R);
		double u = params.M_inf * params.a_inf;
		double Re = rho * u * params.plate_length / params.mu;
		double delta = 5 * params.plate_length / sqrt(Re);

		double lvert = 5 * delta;
		return lvert / size;
	}

	void SupersonicRocketNozzle::InitializeFlowFieldVariables() {

		for (int i = 0; i < imax_; i++) {
			for (int j = 0; j < jmax_; j++) {
				T_.Get(i, j) = flow_parameters_.T_inf;
				P_.Get(i, j) = flow_parameters_.P_inf;
				rho_.Get(i, j) = flow_parameters_.P_inf / flow_parameters_.T_inf / flow_parameters_.R;
				u_.Get(i, j) = 0;
				v_.Get(i, j) = 0;
				M_.Get(i, j) = 0;
				e_.Get(i, j) = flow_parameters_.cv * flow_parameters_.T_inf;

				outside_.Set(i, j, INSIDE);

				if (i < nozzle_exit_width && j <= nozzle_exit_height) {
					outside_.Set(i, j, OUTSIDE);
				}

				if (j < nozzle_exit_height + nozzle_wall_height && j > nozzle_exit_height &&
					i < nozzle_wall_width) {
					outside_.Set(i, j, OUTSIDE);
				}
			}
		}

	}

	void SupersonicRocketNozzle::BoundaryConditions(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e) {

		InletBoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

		NozzleWallBoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

		NozzleExitBoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

		SymmetryBoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

		TopBoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

		OutletBoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

	}

	void SupersonicRocketNozzle::InletBoundaryConditions(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e) {

		for (int j = JMAX - 1; j >= JMAX - freestream_inlet; j--) {
			u.Set(0, j, 0);
			v.Set(0, j, 0);
			P.Set(0, j, params.P_inf);
			T.Set(0, j, params.T_inf);
			rho.Set(0, j, params.P_inf / params.R / params.T_inf);
			e.Set(0, j, params.T_inf * params.cv);

			outside_.Set(0, j, BOUNDARY);
		}
	}

	void SupersonicRocketNozzle::NozzleExitBoundaryConditions(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e) {

		for (int j = 0; j < nozzle_exit_height; j++) {
			u.Set(nozzle_exit_width, j, params.a_inf * params.M_inf);
			v.Set(nozzle_exit_width, j, 0);
			P.Set(nozzle_exit_width, j, P_exit_);
			T.Set(nozzle_exit_width, j, T_exit_);
			rho.Set(nozzle_exit_width, j, P_exit_ / params.R / T_exit_);
			e.Set(nozzle_exit_width, j, T_exit_ * params.cv);

			outside_.Set(nozzle_exit_width, j, BOUNDARY);
		}
	}

	void SupersonicRocketNozzle::NozzleWallBoundaryConditions(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e) {

		double p_ij;
		for (int i = nozzle_exit_width; i < nozzle_wall_width; i++) {
			p_ij = 2 * P.Get(i, nozzle_exit_height - 1) - P.Get(i, nozzle_exit_height - 2);
			T_.Get(i, nozzle_exit_height) = flow_parameters_.T_wall;
			u_.Get(i, nozzle_exit_height) = 0;
			v_.Get(i, nozzle_exit_height) = 0;
			P.Set(i, nozzle_exit_height, p_ij);
			rho_.Get(i, nozzle_exit_height) = p_ij / flow_parameters_.T_wall / flow_parameters_.R;
			e_.Get(i, nozzle_exit_height) = flow_parameters_.cv * flow_parameters_.T_wall;
		
			outside_.Set(i, nozzle_exit_height, BOUNDARY);
		}

		for (int j = nozzle_exit_height; j < nozzle_exit_height + nozzle_wall_height; j++) {
			p_ij = 2 * P.Get(nozzle_wall_width + 1, j) - P.Get(nozzle_wall_width + 2, j);
			T_.Get(nozzle_wall_width, j) = flow_parameters_.T_wall;
			u_.Get(nozzle_wall_width, j) = 0;
			v_.Get(nozzle_wall_width, j) = 0;
			P.Set(nozzle_wall_width, j, p_ij);
			rho_.Get(nozzle_wall_width, j) = p_ij / flow_parameters_.T_wall / flow_parameters_.R;
			e_.Get(nozzle_wall_width, j) = flow_parameters_.cv * flow_parameters_.T_wall;
		
			outside_.Set(nozzle_wall_width, j, BOUNDARY);
		}

		for (int i = 0; i <= nozzle_wall_width; i++) {
			p_ij = 2 * P.Get(i, nozzle_exit_height + nozzle_exit_height + nozzle_wall_height + 1) - 
				P.Get(i, nozzle_exit_height + nozzle_wall_height + 2);
			T_.Get(i, nozzle_exit_height + nozzle_wall_height) = flow_parameters_.T_wall;
			u_.Get(i, nozzle_exit_height + nozzle_wall_height) = 0;
			v_.Get(i, nozzle_exit_height + nozzle_wall_height) = 0;
			P.Set(i, nozzle_exit_height + nozzle_wall_height, p_ij);
			rho_.Get(i, nozzle_exit_height + nozzle_wall_height) = p_ij / flow_parameters_.T_wall / flow_parameters_.R;
			e_.Get(i, nozzle_exit_height + nozzle_wall_height) = flow_parameters_.cv * flow_parameters_.T_wall;
		
			outside_.Set(i, nozzle_exit_height + nozzle_wall_height, BOUNDARY);
		}
	}

	void SupersonicRocketNozzle::SymmetryBoundaryConditions(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e) {

		double p_ij;
		double rho_ij;
		double u_ij;
		double T_ij;
		double K;
		for (int i = nozzle_exit_width + 1; i < IMAX; i++) {
			K = 2; // r(3) / r(2) = 2*delta_y / delta_y
			rho_ij = (K * K * rho.Get(i, 1) - rho.Get(i, 2)) / (K * K - 1);
			u_ij = (K * K * u.Get(i, 1) - u.Get(i, 2)) / (K * K - 1);
			p_ij = 2 * P.Get(i, 1) - P.Get(i, 2);

			// c_p*T + u*u/2 is constant
			//T_ij = (flow_parameters_.cp*T.Get(i - 1, 0) + u.Get(i - 1, 0)* u.Get(i - 1, 0) / 2 - 
			//	u.Get(i, 0)* u.Get(i, 0) / 2) / flow_parameters_.cp;
			//p_ij = T_ij * params.R * rho_ij;
			T_ij = p_ij / params.R / rho_ij;

			u.Set(i, 0, u_ij);
			v.Set(i, 0, 0);
			rho.Set(i, 0, rho_ij);
			P.Set(i, 0, p_ij);
			T.Set(i, 0, T_ij);
			e.Set(i, 0, T_ij * params.cv);

			outside_.Set(i, 0, BOUNDARY);
		}
	}

	void SupersonicRocketNozzle::TopBoundaryConditions(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e) {

		for (int j = 0; j < jmax; j++) {
			u.Set(imax - 1, j, 2 * u.Get(imax - 2, j) - u.Get(imax - 3, j));
			v.Set(imax - 1, j, 2 * v.Get(imax - 2, j) - v.Get(imax - 3, j));
			P.Set(imax - 1, j, 2 * P.Get(imax - 2, j) - P.Get(imax - 3, j));
			T.Set(imax - 1, j, 2 * T.Get(imax - 2, j) - T.Get(imax - 3, j));
			rho.Set(imax - 1, j, 2 * rho.Get(imax - 2, j) - rho.Get(imax - 3, j));
			e.Set(imax - 1, j, 2 * e.Get(imax - 2, j) - e.Get(imax - 3, j));;

			outside_.Set(imax - 1, j, BOUNDARY);
		}
	}

	void SupersonicRocketNozzle::OutletBoundaryConditions(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e) {

		for (int i = 0; i < imax; i++) {
			u.Set(i, jmax - 1, 2 * u.Get(i, jmax - 2) - u.Get(i, jmax - 3));
			v.Set(i, jmax - 1, 2 * v.Get(i, jmax - 2) - v.Get(i, jmax - 3));
			P.Set(i, jmax - 1, 2 * P.Get(i, jmax - 2) - P.Get(i, jmax - 3));
			T.Set(i, jmax - 1, 2 * T.Get(i, jmax - 2) - T.Get(i, jmax - 3));
			rho.Set(i, jmax - 1, 2 * rho.Get(i, jmax - 2) - rho.Get(i, jmax - 3));
			e.Set(i, jmax - 1, 2 * e.Get(i, jmax - 2) - e.Get(i, jmax - 3));

			outside_.Set(i, jmax - 1, BOUNDARY);
		}
	}

}