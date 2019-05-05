#pragma once
#include "supersonic_flow.h"
#include "flow_parameters.h"
#include "math_utils.h"
#include "maccormack.h"


namespace supersonic_rocket_nozzle {
	class SupersonicRocketNozzle : public supersonic_flow::SupersonicFlow
	{
	public:
		SupersonicRocketNozzle();
		~SupersonicRocketNozzle();

		void InitializeFlowFieldVariables();

		void BoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

	private:
		double CalcXStep(const FlowParameters& params, int size);
		double CalcYStep(const FlowParameters& params, int size);

		void InletBoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

		void NozzleExitBoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

		void NozzleWallBoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

		void SymmetryBoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

		void OutletBoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

		void TopBoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

	private:
		double T_exit_;
		double P_exit_;
	};
}

