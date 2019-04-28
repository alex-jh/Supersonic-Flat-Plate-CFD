#pragma once
#include "flow_parameters.h"
#include "math_utils.h"
#include "maccormack.h"

namespace supersonic_plate {
	class SupersonicPlate
	{
	public:
		SupersonicPlate();
		~SupersonicPlate();

		void Run();
		void InitializeFlowFliedVariables();

		void BoundaryConditions(int imax, int jmax, const FlowParameters& params,
			Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
			Array2D<double>& T, Array2D<double>& e);

		bool CheckConvergence(double& diff);

	private:
		double CalcXStep(const FlowParameters& params, int size);
		double CalcYStep(const FlowParameters& params, int size);

	private:
		FlowParameters flow_parameters_;
		MacCormack maccormack_solver_;
		int imax_;
		int jmax_;
		int maxit;
		double deltax_;
		double deltay_;

		Array2D<double> T_;
		Array2D<double> u_;
		Array2D<double> v_;
		Array2D<double> P_;
		Array2D<double> rho_;
		Array2D<double> M_;
		Array2D<double> e_;

		Array2D<double> T_old_;
		Array2D<double> u_old_;
		Array2D<double> v_old_;
		Array2D<double> P_old_;
		Array2D<double> rho_old_;
		Array2D<double> e_old_;

		Array2D<NODE_TYPE> outside_;
	};
}
