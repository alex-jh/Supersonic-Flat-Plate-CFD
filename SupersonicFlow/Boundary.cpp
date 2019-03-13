#include "Boundary.h"
#include "State.h"


Boundary::Boundary() {
}


Boundary::~Boundary()
{
}

NodeState DirichletBoundary::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {
	return last_state[y][x];
}
