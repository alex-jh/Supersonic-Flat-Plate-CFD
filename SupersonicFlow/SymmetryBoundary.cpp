#include "SymmetryBoundary.h"
#include "State.h"
#include "FlowUtils.h"
#include <assert.h>

SymmetryBoundary::SymmetryBoundary()
{
}


SymmetryBoundary::~SymmetryBoundary()
{
}

NodeState SymmetryBoundary::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {

	NodeState state = last_state[y][x];

	state[2] = 0.0;

	double mach = flow_utils::calcMach(last_state[y][x]);

	int dy = 1;
	if (y > last_state.getYSize() / 2) dy = -1;

	if (mach >= 1) {
		double K = abs(last_state.getY(y + 2 * dy, x) - last_state.getY(y, x)) /
			abs(last_state.getY(y + dy, x) - last_state.getY(y, x));

		state[0] = (K * K*last_state[y + dy][x][0] - last_state[y + 2 * dy][x][0]) /
			(K*K - 1);
		state[1] = (K * K*last_state[y + dy][x][1] / last_state[y + dy][x][0] - 
			last_state[y + 2 * dy][x][1] / last_state[y + 2 * dy][x][1]) /
			(K*K - 1) * state[0];
	}
	else {
		state[0] = last_state[y + dy][x][0];
		state[1] = last_state[y + dy][x][1];
	}

	if (cur_state[y][x - 1].done) {

		double T = flow_utils::calcTemperature(cur_state[y][x - 1]);

		double enthalpy = T * flow_utils::getCp(cur_state[y][x - 1]) + pow(cur_state[y][x - 1][1] / cur_state[y][x - 1][0], 2) / 2;
		T = (enthalpy - pow(state[1] / state[0], 2) / 2) / flow_utils::getCp(state);

		state[3] = T * flow_utils::getCv(state) + pow(state[1] / state[0], 2) / 2;
		state[3] *= state[0];
	}

	state.done = true;

	return state;
}
