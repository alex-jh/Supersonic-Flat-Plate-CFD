#pragma once

class NodeState;
class State;
class Solver;

// The order is important
enum class BoundaryType {
	EMPTY = 0,
	WALL_BOUNDARY,
	DIRICHLET_BOUNDARY,
	UNKNOWN_BOUNDARY,
	CAVITY_BOUNDARY,
	CAVITY_BOUNDARY_X,
	CAVITY_BOUNDARY_Y,
	PRIORITY_1_BOUNDARY,
	SYMMETRY_BOUNDARY,
	MACH_NEUMANN_BOUNDARY,
	NUM_BOUNDARIES,
};

class Boundary
{
public:
	Boundary();
	~Boundary();

	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t=1.0) = 0;

	virtual BoundaryType getType() = 0;
};

class DirichletBoundary : public Boundary {
public:
	virtual NodeState calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t = 1.0);

	BoundaryType getType() { return BoundaryType::DIRICHLET_BOUNDARY; }
};
