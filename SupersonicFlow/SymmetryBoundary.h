#pragma once
#include "Boundary.h"

class SymmetryBoundary : public Boundary
{
public:
	SymmetryBoundary();
	~SymmetryBoundary();

	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t=1.0);

	BoundaryType getType() { return BoundaryType::SYMMETRY_BOUNDARY; }
};

