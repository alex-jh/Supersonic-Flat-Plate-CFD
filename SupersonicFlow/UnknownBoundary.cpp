#include "UnknownBoundary.h"
#include "FlowUtils.h"
#include <assert.h>


UnknownBoundary::UnknownBoundary()
{
}


UnknownBoundary::~UnknownBoundary()
{
}

NodeState UnknownBoundary::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {

	double mach = flow_utils::calcMach(last_state[y][x]);

	int dx = 0;
	int dy = 0;
	if (x < last_state.getXSize() / 2) dx = 1;
	else if (x > last_state.getXSize() / 2) dx = -1;
	else if (y < last_state.getYSize() / 2) dy = 1;
	else if (y > last_state.getYSize() / 2) dx = -1;
	else {
		assert(0);
	}

	NodeState state;
	if (mach >= 1.0) {
		state = last_state[y + dy][x + dx] * 3 - last_state[y + 2 * dy][x + 2 * dx] * 3
			+ last_state[y + 3 * dy][x + 3 * dx];
	}
	else {
		state = last_state[y + dy][x + dx];
	}

	state.boundary_ = last_state[y][x].boundary_;

	return state;
}