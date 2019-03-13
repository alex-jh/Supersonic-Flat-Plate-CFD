#pragma once

class State;
class NodeState;

typedef enum IntegrationType_ {
	EULER_BACKWARD = 0,
	EULER_FORWARD,
	CENTRAL_DIFFERENCES
} IntegrationType;

class Solver
{
public:
	Solver();
	~Solver();

	static Solver& GetSolver() {
		static Solver solver;
		return solver;
	}

	double getDeltaTx();

	State solve(State& last_state);
	State solvePredictorStep(State& last_state, int delta, NodeState(Solver::*func)(int, int, State&));
	State solveCorrectorStep(State& last_state, State& predictor_state, int delta, NodeState(Solver::*func)(int, int, State&));
	NodeState solveInternalPredictorStep(int y, int x, double delta_t, double delta, State& last_state, NodeState(Solver::*func)(int, int, State&));
	NodeState solveInternalCorrectorStep(int y, int x, double delta_t, double delta, State& last_state, State& predictor_state, NodeState(Solver::*func)(int, int, State&));

	State solvePredictorStep(State& last_state, int mask = 3, double delta_multiplier = 1.0);
	State solveCorrectorStep(State& last_state, State& predictor_state, int mask = 3, double delta_multiplier = 1.0);

	NodeState solveInternalPredictorStep(int y, int x, double delta_t, State& last_state, int mask = 3);
	NodeState solveInternalCorrectorStep(int y, int x, double delta_t, State& last_state, State& predictor_state, int mask = 3);

	void calcBoundaryNodes(State& last_state, State& state, double delta_tx);

	NodeState F(int y, int x, State& state, IntegrationType integrationType);
	NodeState G(int y, int x, State& state, IntegrationType integrationType);
	NodeState G2(int y, int x, State& state, IntegrationType integrationType);
	NodeState H(int y, int x, State& state, IntegrationType integrationType);
	NodeState J(int y, int x, State& state, IntegrationType integrationType);

	double getDerivativeX(State& state, int x, int y, IntegrationType integrationType, double(*func)(State&, int, int));
	double getDerivativeY(State& state, int x, int y, IntegrationType integrationType, double(*func)(State&, int, int));

	void VerifyState(State& state);

	bool axilsymmetric;
	double delta_t;
};
