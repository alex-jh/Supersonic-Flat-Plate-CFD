#include "CavityBoundaries.h"
#include "State.h"
#include "Solver.h"
#include "FlowUtils.h"

NodeState ConstantVelocityAndPressure::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {
	return last_state[y][x];
}

NodeState ConstantVelocityAndPressureX::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {
	if (cur_state[y][x].size() > 0) {
		return cur_state[y][x];
	}
	
	NodeState state(4);

	double rho = last_state[y][x][0];
	double u = last_state[y][x][1] / last_state[y][x][0];
	double v = last_state[y][x][2] / last_state[y][x][0];

	state[0] = rho;
	state[1] = rho * u;
	state[2] = rho * v;

	bool found = false;
	double P;
	if (x < last_state.getXSize() / 2) {
		for (int i = x + 1; i < x + 5; i++) {
			if (cur_state[y][i].size() > 0) {
				P = flow_utils::calcPressure(cur_state[y][i]);
				found = true;
				break;
			}
		}
	}
	else {
		for (int i = x - 1; i >= x-5; i--) {
			if (cur_state[y][i].size() > 0) {
				P = flow_utils::calcPressure(cur_state[y][i]);
				found = true;
				break;
			}
		}
	}

	if (!found) {
		return NodeState();
	}

	state[3] = flow_utils::calcEnergyWithPressure(state, P);

	return state;
}

NodeState ConstantVelocityAndPressureY::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {
	if (cur_state[y][x].size() > 0) {
		return cur_state[y][x];
	}
	
	NodeState state(4);

	double rho = last_state[y][x][0];
	double u = last_state[y][x][1] / last_state[y][x][0];
	double v = last_state[y][x][2] / last_state[y][x][0];

	state[0] = rho;
	state[1] = rho * u;
	state[2] = rho * v;

	double P = -1;
	bool found = false;
	if (y < last_state.getYSize() / 2) {
		for (int i = y + 1; i < y+5; i++) {
			if (cur_state[i][x].size() > 0) {
				P = flow_utils::calcPressure(cur_state[i][x]);
				found = true;
				break;
			}
		}
	}
	else {
		for (int i = y - 1; i >= y-5; i--) {
			if (cur_state[i][x].size() > 0) {
				P = flow_utils::calcPressure(cur_state[i][x]);
				found = true;
				break;
			}
		}
	}

	if (!found) {
		return NodeState();
	}

	state[3] = flow_utils::calcEnergyWithPressure(state, P);

	return state;
}