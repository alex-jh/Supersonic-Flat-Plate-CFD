#pragma once
#include "Boundary.h"

class WallBoundary : public Boundary
{
public:
	WallBoundary();
	~WallBoundary();

	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t = 1.0);

	BoundaryType getType() { return BoundaryType::WALL_BOUNDARY; }
};

