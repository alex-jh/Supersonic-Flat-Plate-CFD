#include "Solver.h"
#include "State.h"
#include "FlowUtils.h"
#include "Boundary.h"
#include "MachNeumannBoundary.h"
#include <assert.h>
#include <assert.h>

bool pressure_in_conv_energy = true;
bool debug = false;

double Solver::getDeltaTx() {
	return delta_t;
}

Solver::Solver() {
	axilsymmetric = true;
	delta_t = 0.1;
}


Solver::~Solver()
{
}

State Solver::solve(State& last_state) {
	bool split = true;
	if (split) {
		State predictor_state = solvePredictorStep(last_state, 2, 0.5);
		State state = solveCorrectorStep(last_state, predictor_state, 2, 0.5);

		predictor_state = solvePredictorStep(state, 1);
		state = solveCorrectorStep(state, predictor_state, 1);

		predictor_state = solvePredictorStep(state, 2, 0.5);
		state = solveCorrectorStep(state, predictor_state, 2, 0.5);
	
		return state;
	}
	else {
		State predictor_state = solvePredictorStep(last_state);
		State state = solveCorrectorStep(last_state, predictor_state);
		
		return state;
	}
}

State Solver::solvePredictorStep(State& last_state, int delta, NodeState(Solver::*func)(int, int, State&)) {

	double delta_t = getDeltaTx() / 2;

	State predictor_state(State::createStateFrom(last_state));

	int deltax = 1;
	int deltay = 0;

	if (delta == -1) {
		deltax = 0;
		deltay = 1;
	}

	// The internal points.
	// Boundary points with Neumann type, we have found the nodes outside
	// the boundary, and boundary points are calculated the usual way.
	for (int y = 1; y < last_state.getYSize() - 1; y++) {
		for (int x = 1; x < last_state.getXSize() - 1; x++) {
			if (!last_state.isWall(y, x) && (last_state[y][x].boundary_ == NULL ||
					last_state[y][x].boundaryType_ == BoundaryType::MACH_NEUMANN_BOUNDARY)) {
				double delta_s = abs(last_state.getX(y, x) - last_state.getX(y + deltay, x + deltax));
				predictor_state[y][x] = solveInternalPredictorStep(y, x, delta_t, delta_s, last_state, func);
			}
		}
	}

	// Calculate boundary nodes.
	calcBoundaryNodes(last_state, predictor_state, delta_t);

	return predictor_state;
}

State Solver::solvePredictorStep(State& last_state, int mask, double delta_multiplier) {

	double delta_tx = getDeltaTx() * delta_multiplier;

	State predictor_state(State::createStateFrom(last_state));

	// The internal points.
	// Boundary points with Neumann type, we have found the nodes outside
	// the boundary, and boundary points are calculated the usual way.
	for (int y = 1; y < last_state.getYSize() - 1; y++) {
		for (int x = 1; x < last_state.getXSize() - 1; x++) {
			if (!last_state.isWall(y, x) && (last_state[y][x].boundary_ == NULL ||
				last_state[y][x].boundaryType_ == BoundaryType::MACH_NEUMANN_BOUNDARY)) {
				predictor_state[y][x] = solveInternalPredictorStep(y, x, delta_tx, last_state, mask);
				predictor_state[y][x].boundary_ = last_state[y][x].boundary_;
			}
		}
	}

	// Calculate boundary nodes.
	calcBoundaryNodes(last_state, predictor_state, delta_tx);

	return predictor_state;
}

State Solver::solveCorrectorStep(State& last_state, State& predictor_state, int delta, NodeState(Solver::*func)(int, int, State&)) {

	double delta_t = getDeltaTx() / 2;

	State corrector_state = last_state;

	int deltax = 1;
	int deltay = 0;

	if (delta == -1) {
		deltax = 0;
		deltay = 1;
		delta_t = getDeltaTx();
	}

	// The internal points.
	// Boundary points with Neumann type, we have found the nodes outside
	// the boundary, and boundary points are calculated the usual way.
	for (int y = 0; y < last_state.getYSize(); y++) {
		for (int x = 0; x < last_state.getXSize(); x++) {
			if (!last_state.isWall(y, x) && (last_state[y][x].boundary_ == NULL ||
				last_state[y][x].boundaryType_ == BoundaryType::MACH_NEUMANN_BOUNDARY)) {
				double delta_s = abs(last_state.getX(y, x) - last_state.getX(y + deltay, x + deltax));
				corrector_state[y][x] = solveInternalCorrectorStep(y, x, delta_t, delta_s, last_state, predictor_state, func);
			}
		}
	}

	// Calculate boundary nodes.
	calcBoundaryNodes(last_state, corrector_state, delta_t);

	return corrector_state;
}

State Solver::solveCorrectorStep(State& last_state, State& predictor_state, int mask, double delta_multiplier) {

	double delta_tx = getDeltaTx() * delta_multiplier;

	State corrector_state = last_state;

	// the internal points.
	// Boundary points with Neumann type, we have found the nodes outside
	// the boundary, and boundary points are calculated the usual way.
	for (int y = 0; y < last_state.getYSize(); y++) {
		for (int x = 0; x < last_state.getXSize(); x++) {
			if (!last_state.isWall(y, x) && (last_state[y][x].boundary_ == NULL ||
				last_state[y][x].boundaryType_ == BoundaryType::MACH_NEUMANN_BOUNDARY)) {
				corrector_state[y][x] = solveInternalCorrectorStep(y, x, delta_tx, last_state, predictor_state, mask);
				corrector_state[y][x].boundary_ = last_state[y][x].boundary_;
			}
		}
	}

	calcBoundaryNodes(predictor_state, corrector_state, delta_tx);

	return corrector_state;
}

void Solver::calcBoundaryNodes(State& last_state, State& state, double delta_tx) {
	
	std::vector<std::pair<NodeState*, std::pair<int, int> > > boundaryNodes;
	boundaryNodes.reserve(5 * (state.getXSize() + state.getYSize()));

	for (int y = 0; y < last_state.getYSize(); y++) {
		for (int x = 0; x < last_state.getXSize(); x++) {
			if (last_state[y][x].boundary_ != NULL) {
				boundaryNodes.push_back(std::make_pair(
					&last_state[y][x], std::make_pair(y, x)));
			}
		}
	}

	std::sort(boundaryNodes.begin(), boundaryNodes.end(),
		[](const std::pair<NodeState*, std::pair<int, int>>& a,
			const std::pair<NodeState*, std::pair<int, int>>& b) -> bool
	{
		return a.first->boundaryType_ < b.first->boundaryType_;
	});

	for (BoundaryType boundaryType = BoundaryType::DIRICHLET_BOUNDARY;
		boundaryType < BoundaryType::MACH_NEUMANN_BOUNDARY; 
		boundaryType = (BoundaryType)((int)boundaryType + 1)) {
		for (auto& nodes : boundaryNodes) {
			for (auto* boundary : *nodes.first->boundary_) {
				if (boundary->getType() == boundaryType) {
					int y = nodes.second.first;
					int x = nodes.second.second;
					state[y][x] =
						boundary->calcState(x, y, last_state, *this, state, delta_tx);
					state[y][x].boundary_ = last_state[y][x].boundary_;
				}
			}
		}
	}

	//for (int y = 0; y < last_state.getYSize(); y++) {
	//	for (int x = 0; x < last_state.getXSize(); x++) {
	//		if (last_state[y][x].boundary_ != NULL &&
	//				last_state[y][x].boundaryType_ <= BoundaryType::PRIORITY_1_BOUNDARY) {
	//			for (auto* boundary : *last_state[y][x].boundary_) {
	//				state[y][x] =
	//					boundary->calcState(x, y, last_state, *this, state, delta_tx);
	//				state[y][x].boundary_ = last_state[y][x].boundary_;
	//			}
	//		}
	//	}
	//}

	//for (int y = 0; y < last_state.getYSize(); y++) {
	//	for (int x = 0; x < last_state.getXSize(); x++) {
	//		if (last_state[y][x].boundary_ != NULL &&
	//			last_state[y][x].boundaryType_ < BoundaryType::MACH_NEUMANN_BOUNDARY) {
	//			for (auto* boundary : *last_state[y][x].boundary_) {
	//				state[y][x] =
	//					boundary->calcState(x, y, last_state, *this, state, delta_tx);
	//				state[y][x].boundary_ = last_state[y][x].boundary_;
	//			}
	//		}
	//	}
	//}

	if (debug) {
		VerifyState(state);
	}

	MachNeumannBoundary::calcExternalBoundary(last_state);
}

NodeState Solver::solveInternalPredictorStep(int y, int x, double delta_t, double delta, State& last_state, NodeState(Solver::*func)(int, int, State&)) {

	NodeState predictor_state = last_state[y][x];

	// If used correct the delta it should be either (y, x - 1, last_state) or (y - 1, x, last_state). Pass deltax and deltay tothe function.
	predictor_state -= delta_t / delta * ((this->*func)(y, x, last_state) - (this->*func)(y, x - 1, last_state)) +
		delta_t / last_state.getY(y, x) * H(y, x, last_state, EULER_FORWARD);

	return predictor_state;
}

NodeState Solver::solveInternalPredictorStep(int y, int x, double delta_t, State& last_state, int mask) {

	NodeState predictor_state = last_state[y][x];

	double delta_x = abs(last_state.getX(y, x) - last_state.getX(y, x - 1));
	double delta_y = abs(last_state.getY(y, x) - last_state.getY(y - 1, x));

	//NodeState state_x = -delta_t / delta_x * (F(y, x, last_state) - F(y, x - 1, last_state));
	//NodeState state_y = -delta_t / delta_y * (last_state.getY(y, x)*G(y, x, last_state) - last_state.getY(y - 1, x)*G(y - 1, x, last_state)) /
	//	last_state.getY(y, x);
	//NodeState state_h = delta_t * (0.5*H(y, x, last_state) / last_state.getY(y, x) + 0.5*H(y - 1, x, last_state) / last_state.getY(y, x));

	NodeState state_x(4);
	if (mask & 2) {
		state_x = -delta_t / delta_x * (F(y, x, last_state, EULER_FORWARD) - F(y, x - 1, last_state, EULER_FORWARD));
	}
	NodeState state_y(4);
	if (mask & 1) {
		if (axilsymmetric) {
			state_y = -delta_t / delta_y * (last_state.getY(y, x)*G2(y, x, last_state, EULER_FORWARD) - last_state.getY(y - 1, x)*G2(y - 1, x, last_state, EULER_FORWARD)) /
				last_state.getY(y, x) - delta_t / delta_y * (J(y, x, last_state, EULER_FORWARD) - J(y - 1, x, last_state, EULER_FORWARD));
		}
		else {
			state_y = -delta_t / delta_y * (G2(y, x, last_state, EULER_FORWARD) + J(y, x, last_state, EULER_FORWARD) 
				- G2(y - 1, x, last_state, EULER_FORWARD) - J(y - 1, x, last_state, EULER_FORWARD));
		}
	}
	NodeState state_h(4);

	predictor_state += state_x + state_y + state_h;

	double T = flow_utils::calcTemperature(predictor_state);
	double T1 = flow_utils::calcTemperature(last_state[y][x-1]);
	double T2 = flow_utils::calcTemperature(last_state[y][x]);

	return predictor_state;
}

NodeState Solver::solveInternalCorrectorStep(int y, int x, double delta_t, double delta, State& last_state, State& predictor_state, NodeState(Solver::*func)(int, int, State&)) {

	NodeState corrector_state = last_state[y][x] + predictor_state[y][x];

	corrector_state = - delta_t / delta * ((this->*func)(y, x + 1, predictor_state) - (this->*func)(y, x, predictor_state)) +
		delta_t / last_state.getY(y, x) * H(y, x, last_state, EULER_BACKWARD);
	corrector_state = corrector_state * 0.5;

	return corrector_state;
}

NodeState Solver::solveInternalCorrectorStep(int y, int x, double delta_t, State& last_state, State& predictor_state, int mask) {

	NodeState corrector_state = last_state[y][x] + predictor_state[y][x];

	double delta_x = abs(last_state.getX(y, x + 1) - last_state.getX(y, x));
	double delta_y = abs(last_state.getY(y + 1, x) - last_state.getY(y, x));

	//NodeState state_x = -delta_t / delta_x * (F(y, x + 1, predictor_state) - F(y, x, predictor_state));
	//NodeState state_y = -delta_t / delta_y * (last_state.getY(y + 1, x)*G(y + 1, x, predictor_state) - last_state.getY(y, x)*G(y, x, predictor_state)) /
	//	last_state.getY(y, x);
	//NodeState state_h = delta_t * (0.5*H(y + 1, x, predictor_state) / last_state.getY(y, x) + 0.5*H(y, x, predictor_state) / last_state.getY(y, x));

	NodeState state_x(4);
	if (mask & 2) {
		state_x = -delta_t / delta_x * (F(y, x + 1, predictor_state, EULER_BACKWARD) - F(y, x, predictor_state, EULER_BACKWARD));
	}
	NodeState state_y(4);
	if (mask & 1) {
		if (axilsymmetric) {
			state_y = -delta_t / delta_y * (last_state.getY(y + 1, x)*G2(y + 1, x, predictor_state, EULER_BACKWARD) - last_state.getY(y, x)*G2(y, x, predictor_state, EULER_BACKWARD)) /
				last_state.getY(y, x) - delta_t / delta_y * (J(y + 1, x, predictor_state, EULER_BACKWARD) - J(y, x, predictor_state, EULER_BACKWARD));
		}
		else {
			state_y = -delta_t / delta_y * (G2(y + 1, x, predictor_state, EULER_BACKWARD) + J(y + 1, x, predictor_state, EULER_BACKWARD)
				- G2(y, x, predictor_state, EULER_BACKWARD) - J(y, x, predictor_state, EULER_BACKWARD));
		}
	}
	NodeState state_h(4);

	corrector_state += state_x + state_y + state_h;
	corrector_state = corrector_state * 0.5;

	return corrector_state;
}

NodeState Solver::F(int y, int x, State& state, IntegrationType integrationType) {
	NodeState f(4);

	f[0] = state[y][x][1];

	double u = state[y][x][1] / state[y][x][0];
	double v = state[y][x][2] / state[y][x][0];

	double r = state.getY(y, x);

	double u_x = getDerivativeX(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return state[y][x][1] / state[y][x][0];
	});
	double u_r = getDerivativeY(state, x, y, CENTRAL_DIFFERENCES, [](State& state, int y, int x) -> double
	{
		return state[y][x][1] / state[y][x][0];
	});

	double v_x = getDerivativeX(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return state[y][x][2] / state[y][x][0];
	});
	double v_r = getDerivativeY(state, x, y, CENTRAL_DIFFERENCES, [](State& state, int y, int x) -> double
	{
		return state[y][x][2] / state[y][x][0];
	});

	double div_v = u_x + v_r + axilsymmetric * (v / r);

	double P = flow_utils::calcPressure(state[y][x]);
	double mu = flow_utils::calcMu(state[y][x]);
	double eps = flow_utils::calcEps(state[y][x]);
	double lambda = -2.0 / 3.0*(mu + eps);

	double sigma_xx = lambda * div_v + 2 * (mu + eps)*u_x;

	f[1] = u * u*state[y][x][0] - sigma_xx + P;

	double tau_xr = (mu + eps)*(u_r + v_x);

	f[2] = u * v*state[y][x][0] - tau_xr;

	double T = flow_utils::calcTemperature(state[y][x]);
	double Pr = flow_utils::getPr();
	double Pr_t = flow_utils::getPr_t();
	double T_x = getDerivativeX(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return flow_utils::calcTemperature(state[y][x]);
	});
	double q_x = -flow_utils::getCp(T)*(mu / Pr + eps / Pr_t)*T_x;

	if (pressure_in_conv_energy) {
		sigma_xx -= P;
	}

	f[3] = u * state[y][x][3] + q_x - u * sigma_xx - v * tau_xr;

	return f;
}

NodeState Solver::G(int y, int x, State& state, IntegrationType integrationType) {
	NodeState g(4);

	g[0] = state[y][x][2];

	double u = state[y][x][1] / state[y][x][0];
	double v = state[y][x][2] / state[y][x][0];

	double r = state.getY(y, x);

	double u_x = getDerivativeX(state, x, y, CENTRAL_DIFFERENCES, [](State& state, int y, int x) -> double
	{
		return state[y][x][1] / state[y][x][0];
	});
	double u_r = getDerivativeY(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return state[y][x][1] / state[y][x][0];
	});

	double v_x = getDerivativeX(state, x, y, CENTRAL_DIFFERENCES, [](State& state, int y, int x) -> double
	{
		return state[y][x][2] / state[y][x][0];
	});
	double v_r = getDerivativeY(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return state[y][x][2] / state[y][x][0];
	});

	double div_v = u_x + v_r + v / r;

	double P = flow_utils::calcPressure(state[y][x]);
	double mu = flow_utils::calcMu(state[y][x]);
	double eps = flow_utils::calcEps(state[y][x]);
	double lambda = -2.0 / 3.0*(mu + eps);

	double tau_xr = (mu + eps)*(u_r + v_x);

	g[1] = u * v*state[y][x][0] - tau_xr;

	double sigma_rr = -P + lambda * div_v + 2 * (mu + eps)*v_r;

	g[2] = v * v*state[y][x][0] - sigma_rr;

	double T = flow_utils::calcTemperature(state[y][x]);
	double Pr = flow_utils::getPr();
	double Pr_t = flow_utils::getPr_t();
	double T_r = getDerivativeY(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return flow_utils::calcTemperature(state[y][x]);
	});

	double q_r = -flow_utils::getCp(T)*(mu / Pr + eps / Pr_t)*T_r;

	g[3] = v * state[y][x][3] + q_r - v * sigma_rr - u * tau_xr;

	return g;
}

NodeState Solver::G2(int y, int x, State& state, IntegrationType integrationType) {
	NodeState g(4);

	g[0] = state[y][x][2];

	double u = state[y][x][1] / state[y][x][0];
	double v = state[y][x][2] / state[y][x][0];

	double r = state.getY(y, x);

	double u_x = getDerivativeX(state, x, y, CENTRAL_DIFFERENCES, [](State& state, int y, int x) -> double
	{
		return state[y][x][1] / state[y][x][0];
	});
	double u_r = getDerivativeY(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return state[y][x][1] / state[y][x][0];
	});

	double v_x = getDerivativeX(state, x, y, CENTRAL_DIFFERENCES, [](State& state, int y, int x) -> double
	{
		return state[y][x][2] / state[y][x][0];
	});
	double v_r = getDerivativeY(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return state[y][x][2] / state[y][x][0];
	});

	double div_v = u_x + v_r + axilsymmetric * (v / r);

	double P = flow_utils::calcPressure(state[y][x]);
	double mu = flow_utils::calcMu(state[y][x]);
	double eps = flow_utils::calcEps(state[y][x]);
	double lambda = -2.0 / 3.0*(mu + eps);

	double tau_xr = (mu + eps)*(u_r + v_x);

	g[1] = u * v*state[y][x][0] - tau_xr;

	double sigma_rr = lambda * div_v + 2 * (mu + eps)*v_r;

	g[2] = v * v*state[y][x][0] - sigma_rr;

	double T = flow_utils::calcTemperature(state[y][x]);
	double Pr = flow_utils::getPr();
	double Pr_t = flow_utils::getPr_t();
	double T_r = getDerivativeY(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return flow_utils::calcTemperature(state[y][x]);
	});

	double q_r = -flow_utils::getCp(T)*(mu / Pr + eps / Pr_t)*T_r;

	g[3] = v * state[y][x][3] + q_r - v * (sigma_rr) - u * tau_xr;

	return g;
}

NodeState Solver::H(int y, int x, State& state, IntegrationType integrationType) {
	NodeState h(4);

	double v = state[y][x][2] / state[y][x][0];

	double r = state.getY(y, x);

	double u_x = getDerivativeX(state, x, y, CENTRAL_DIFFERENCES, [](State& state, int y, int x) -> double
	{
		return state[y][x][1] / state[y][x][0];
	});
	double v_r = getDerivativeY(state, x, y, integrationType, [](State& state, int y, int x) -> double
	{
		return state[y][x][2] / state[y][x][0];
	});

	double div_v = u_x + v_r + v / r;

	double P = flow_utils::calcPressure(state[y][x]);
	double mu = flow_utils::calcMu(state[y][x]);
	double eps = flow_utils::calcEps(state[y][x]);
	double lambda = -2.0 / 3.0*(mu + eps);

	double sigma_h = -P + lambda * div_v + 2 * (mu + eps)*v / r;
	h[2] = -sigma_h;

	return h;
}

NodeState Solver::J(int y, int x, State& state, IntegrationType integrationType) {
	NodeState h(4);


	double P = flow_utils::calcPressure(state[y][x]);
	h[2] = P;

	double v = state[y][x][2] / state[y][x][0];

	if (pressure_in_conv_energy) {
		h[3] = P*v;
	}

	return h;
}

double Solver::getDerivativeX(State& state, int x, int y, IntegrationType integrationType, double(*func)(State&, int, int)) {
	switch (integrationType)
	{
	case EULER_BACKWARD:
	{
		assert(x - 1 >= 0 && state[y][x - 1].vals.size() != 0);
		double delta_x = abs(state.getX(y, x) - state.getX(y, x - 1));
		return (func(state, y, x) - func(state, y, x - 1)) / delta_x;
	}
	case EULER_FORWARD:
	{
		assert(x + 1 < state.getXSize() && state[y][x + 1].vals.size() != 0);
		double delta_x = abs(state.getX(y, x + 1) - state.getX(y, x));
		return (func(state, y, x + 1) - func(state, y, x)) / delta_x;
	}
	case CENTRAL_DIFFERENCES:
	{
		double delta_x = abs(state.getX(y, x + 1) - state.getX(y, x - 1));
		return (func(state, y, x + 1) - func(state, y, x - 1)) / delta_x;
	}
	default:
		return 0.0;
	}
}

double Solver::getDerivativeY(State& state, int x, int y, IntegrationType integrationType, double(*func)(State&, int, int)) {
	switch (integrationType)
	{
	case EULER_BACKWARD:
	{
		assert(y - 1 >= 0 && state[y - 1][x].vals.size() != 0);
		double delta_y = abs(state.getY(y, x) - state.getY(y - 1, x));
		return (func(state, y, x) - func(state, y - 1, x)) / delta_y;
	}
	case EULER_FORWARD:
	{
		assert(y + 1 < state.getYSize() && state[y + 1][x].vals.size() != 0);
		double delta_y = abs(state.getY(y + 1, x) - state.getY(y, x));
		return (func(state, y + 1, x) - func(state, y, x)) / delta_y;
	}
	case CENTRAL_DIFFERENCES:
	{
		double delta_y = abs(state.getY(y + 1, x) - state.getY(y - 1, x));
		return (func(state, y + 1, x) - func(state, y - 1, x)) / delta_y;
	}
	default:
		return 0.0;
	}
}

void Solver::VerifyState(State& state) {
	for (int y = 1; y < state.getYSize() - 1; y++) {
		for (int x = 1; x < state.getXSize() - 1; x++) {
			assert(state[y][x].size() > 0);
			assert(flow_utils::calcPressure(state[y][x]) > 0.0);
		}
	}
}