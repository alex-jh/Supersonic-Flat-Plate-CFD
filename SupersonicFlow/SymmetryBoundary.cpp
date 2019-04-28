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
	if (last_state[y + dy][x].boundary_ != NULL) dy++;
	if (y > last_state.getYSize() / 2) dy = -dy;

	if (mach >= 1) {
		double K = abs(cur_state.getY(y + 2 * dy, x) - cur_state.getY(y, x)) /
			abs(cur_state.getY(y + dy, x) - cur_state.getY(y, x));

		state[0] = (K * K*cur_state[y + dy][x][0] - cur_state[y + 2 * dy][x][0]) /
			(K*K - 1);
		state[1] = (K * K*cur_state[y + dy][x][1] / cur_state[y + dy][x][0] -
			cur_state[y + 2 * dy][x][1] / cur_state[y + 2 * dy][x][1]) /
			(K*K - 1) * state[0];
	}
	else {
		state[0] = cur_state[y + dy][x][0];
		state[1] = cur_state[y + dy][x][1];
	}

	if (cur_state[y][x - 1].done) {

		double T = flow_utils::calcTemperature(cur_state[y + dy][x]);

		double enthalpy = T * flow_utils::getCp(cur_state[y + dy][x]) + (
			pow(cur_state[y + dy][x][1] / cur_state[y + dy][x][0], 2) + pow(cur_state[y + dy][x][2] / cur_state[y + dy][x][0], 2)) / 2;
		T = (enthalpy - pow(state[1] / state[0], 2) / 2) / flow_utils::getCp(state);

		state[3] = T * flow_utils::getCv(state) + pow(state[1] / state[0], 2) / 2;
		state[3] *= state[0];
	}

	state.done = true;

	return state;
}
