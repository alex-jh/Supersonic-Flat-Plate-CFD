#pragma once
#include "Boundary.h"
#include <algorithm>

class Solver;

// Assumes the boundary condition is horizontal.
class MachNeumannBoundary : public Boundary {
public:
	MachNeumannBoundary();
	~MachNeumannBoundary();

	NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t=1.0);
	void calcStateValue(int x, int y, State& last_state, double& x0, double& y0, NodeState& state0);
	std::pair<NodeState, double> calcOutsideStateNodes(int x, int y, State& last_state);
	static void calcExternalBoundary(State& last_state);

	BoundaryType getType() { return BoundaryType::MACH_NEUMANN_BOUNDARY; }
};

