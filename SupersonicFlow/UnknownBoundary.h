#pragma once
#include "Boundary.h"

class UnknownBoundary : public Boundary
{
public:
	UnknownBoundary();
	~UnknownBoundary();

	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t = 1.0);

	BoundaryType getType() { return BoundaryType::UNKNOWN_BOUNDARY; }
};

