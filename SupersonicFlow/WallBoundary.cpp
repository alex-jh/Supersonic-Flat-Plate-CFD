#include "WallBoundary.h"
#include "State.h"
#include "FlowUtils.h"
#include "Constants.h"
#include <assert.h>


WallBoundary::WallBoundary()
{
}


WallBoundary::~WallBoundary()
{
}

NodeState WallBoundary::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {

	NodeState state = last_state[y][x];

	state[2] = 0.0;
	state[1] = 0.0;

	NodeState nxtState = state;
	int dx = 0;
	int dy = 0;

	if (last_state.isWall(y - 2, x)) {
		dy = 1;
	} 
	else if (last_state.isWall(y + 2, x)) {
		dy = -1;
	}
	else if (last_state.isWall(y, x - 2)) {
		dx = 1;
	}
	else if (last_state.isWall(y, x + 2)) {
		dx = -1;
	}
	else if (last_state.isWall(y + 2, x + 2)) {
		dx = -1;
		dy = -1;
	}
	else if (last_state.isWall(y - 2, x - 2)) {
		dx = +1;
		dy = +1;
	}
	else if (last_state.isWall(y - 2, x + 2)) {
		dx = -1;
		dy = 1;
	}
	else if (last_state.isWall(y + 2, x - 2)) {
		dx = 1;
		dy = +1;
	}
	else {
		assert(0 && "Wrong wall condition");
	}

	double P = flow_utils::calcPressure(last_state[y + dy][x + dx]);

	state[0] = flow_utils::calcDensityWithPressure(state, P, constants::T_wall);
	state[3] = flow_utils::calcEnergy(P, constants::T_wall, 0.0, 0.0);

	// nxtState[3] = flow_utils::calcEnergyWithPressure(nxtState, P, constants::T_wall);
	// nxtState[0] = flow_utils::calcDensityWithPressure(nxtState, P, constants::T_wall);

	return state;
}
