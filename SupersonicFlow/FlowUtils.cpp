#include "FlowUtils.h"

namespace flow_utils {

	// Update with:
	// https://www.ohio.edu/mechanical/thermo/property_tables/air/air_cp_cv.html
	double getCv(NodeState& state) {
		return 0.72 * FlowConstants::GetFlowConstants().T_c /
			FlowConstants::GetFlowConstants().u_c / FlowConstants::GetFlowConstants().u_c;
	}

	double getCv(double T) {
		return 0.72 * FlowConstants::GetFlowConstants().T_c /
			FlowConstants::GetFlowConstants().u_c / FlowConstants::GetFlowConstants().u_c;
	}

	double getCp(NodeState& state) {
		return FlowConstants::GetFlowConstants().Cp * FlowConstants::GetFlowConstants().T_c /
			FlowConstants::GetFlowConstants().u_c / FlowConstants::GetFlowConstants().u_c;
	}

	double getCp(double T) {
		return FlowConstants::GetFlowConstants().Cp * FlowConstants::GetFlowConstants().T_c /
			FlowConstants::GetFlowConstants().u_c / FlowConstants::GetFlowConstants().u_c;
	}

	double getLambda() {
		return 1.4;
	}

	double getR() {
		return 287.058 * FlowConstants::GetFlowConstants().T_c / 
			FlowConstants::GetFlowConstants().u_c / FlowConstants::GetFlowConstants().u_c;
	}

	double getEnergy(NodeState& state) {
		return state[3] / state[0];
	}

	double getVx(NodeState& state) {
		return state[1] / state[0];
	}

	double getVy(NodeState& state) {
		return state[2] / state[0];
	}

	double calcTemperature(NodeState& state) {

		double vx = getVx(state);
		double vy = getVy(state);
		double vel2 = (vx * vx + vy * vy);

		return (getEnergy(state) - vel2 / 2) / getCv(state);

	}

	double calcTemperature(double P, double rho) {

		return P / getR() / rho;

	}

	double calcPressure(NodeState& state) {
		double rho = state[0];
		double T = calcTemperature(state);

		return getR() * T * rho;
	}

	double calcEnergyWithPressure(NodeState& state, double P) {
		double vx = getVx(state);
		double vy = getVy(state);
		double vel2 = (vx * vx + vy * vy);

		double T = P / getR() / state[0];

		return getCv(T)*T + vel2 / 2;
	}

	double calcEnergy(double P, double T, double u, double v) {
		double vx = u;
		double vy = v;
		double vel2 = (vx * vx + vy * vy);

		return getCv(T)*T + vel2 / 2;
	}

	double calcDensityWithPressure(NodeState& state, double P, double T) {
		return P / (getR()*T);
	}

	double calcDensity(double P, double T) {
		return P / (getR()*T);
	}

	double calcMach(NodeState& state) {

		double vx = getVx(state);
		double vy = getVy(state);
		double vel = sqrt(vx*vx + vy * vy);

		double T = calcTemperature(state);

		if (T < 0.0) return 0.0;

		double a = getLambda()*getR()*T;

		if (a > 0.0) {
			return vel / sqrt(a);
		}
		else {
			return 0.0;
		}
	}

	double calcMu(NodeState& state) {
		double T = std::max(calcTemperature(state), 1.0e-8);

		double mu0 = 1.458e-6 / FlowConstants::GetFlowConstants().u_c / 
			FlowConstants::GetFlowConstants().P_c * sqrt(FlowConstants::GetFlowConstants().T_c);

		return mu0 *T*sqrt(T) / (T + 110.4);
	}

	double calcEps(NodeState& state) {
		return 0;
	}

	double getPr() {
		return 0.72;
	}

	double getPr_t() {
		return 0.90;
	}
}