#pragma once
#include "Boundary.h"

class ConstantVelocityAndPressure : public Boundary
{
public:

	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t);

	BoundaryType getType() { return BoundaryType::CAVITY_BOUNDARY; }
};

class ConstantVelocityAndPressureX : public Boundary
{
public:

	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t);

	BoundaryType getType() { return BoundaryType::CAVITY_BOUNDARY_X; }
};

class ConstantVelocityAndPressureY : public Boundary
{
public:

	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t);

	BoundaryType getType() { return BoundaryType::CAVITY_BOUNDARY_Y; }
};