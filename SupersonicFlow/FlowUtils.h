#pragma once
#include "State.h"

namespace flow_utils {

	class FlowConstants
	{
	public:
		FlowConstants() {
			P_c = 1;
			T_c = 1;
			u_c = 1;
			r_c = 1;
			Cp = 1.005;
		}

		inline static FlowConstants& GetFlowConstants() {
			static FlowConstants object;
			return object;
		}

		double P_c;
		double T_c;
		double u_c;
		double r_c;
		double Cp;
	};

	// Update with:
	// https://www.ohio.edu/mechanical/thermo/property_tables/air/air_cp_cv.html
	double getCv(NodeState& state);

	double getCv(double T);

	double getCp(NodeState& state);

	double getCp(double T);

	double getLambda();

	double getR();

	double getEnergy(NodeState& state);

	double getVx(NodeState& state);

	double getVy(NodeState& state);

	double calcTemperature(NodeState& state);

	double calcTemperature(double P, double rho);

	double calcPressure(NodeState& state);

	double calcEnergyWithPressure(NodeState& state, double P);

	double calcEnergy(double P, double T, double u, double v);

	double calcDensityWithPressure(NodeState& state, double P, double T);

	double calcDensity(double P, double T);

	double calcMach(NodeState& state);

	double calcMu(NodeState& state);

	double calcEps(NodeState& state);

	double getPr();

	double getPr_t();
}