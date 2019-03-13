#include "MachNeumannBoundary.h"
#include <iostream>
#include "FlowUtils.h"
#include "MathUtils.h"
#include <assert.h>


MachNeumannBoundary::MachNeumannBoundary() {
}

MachNeumannBoundary::~MachNeumannBoundary() {
}

NodeState MachNeumannBoundary::calcState(int x, int y, State& last_state, Solver& solver, State& cur_state, double delta_t) {
	std::pair<NodeState, double> pState = calcOutsideStateNodes(x, y, last_state);
	return pState.first;
}

std::pair<NodeState, double> MachNeumannBoundary::calcOutsideStateNodes(int x, int y, State& last_state) {

	if (x == 0 || y == 0) {
		std::cout << "Error: Neumann boundary conditions need a padding" << std::endl;
	}

	std::vector<NodeState> vStates;
	std::vector<double> vx;
	std::vector<double> vy;

	int xx[] = { x - 1, x, x + 1 };
	for (int i = 0; i < 3; i++) {
		if (xx[i] != 0 && xx[i] != last_state.getXSize() - 1) {
			vx.push_back(0.0);
			vy.push_back(0.0);
			vStates.push_back(NodeState());
			calcStateValue(xx[i], y, last_state, vx.back(), vy.back(), vStates.back());
		}
	}
	
	NodeState state0;
	double y_val;
	if (vStates.size() == 2) {
		state0 = math_utils::interpolateState(vStates[0], vStates[1], abs(vx[0] - vx[1]), vx[1] - last_state.getX(y, x));
		y_val = math_utils::interpolate(vy[0], vy[1], abs(vx[0] - vx[1]), vx[1] - last_state.getX(y, x));
	} 
	else if (vStates.size() == 3) {
		state0 = math_utils::parabolicInterpolateState(vStates, vx, last_state.getX(y, x));
		y_val = math_utils::parabolicInterpolate(vy, vx, last_state.getX(y, x));
	}
	else {
		assert(0 && "Can't find enough information to calculate boundary point");
		std::cout << "Error: Can't find enough information to calculate boundary point" << std::endl;
	}

	return std::make_pair(state0, y_val);

}

void MachNeumannBoundary::calcStateValue(int x, int y, State& last_state, double& x0, double& y0, NodeState& state0) {

	// Mach(u_ij)
	double mach = flow_utils::calcMach(last_state[y][x]);

	int sign = 1;

	if (y != 1) sign = -1;

	// Mach(u_i-1j)
	NodeState pState1 = last_state.getStateWithBoundary(y + sign, x);
	double mach_1 = flow_utils::calcMach(pState1);
	// Mach(u_i-1j-1)
	NodeState pState0 = last_state.getStateWithBoundary(y + sign, x - 1);
	double mach_0 = flow_utils::calcMach(pState0);
	// Mach(u_i-1j+1)
	NodeState pState2 = last_state.getStateWithBoundary(y + sign, x + 1);
	double mach_2 = flow_utils::calcMach(pState2);

	if ((mach >= mach_1 - math_utils::eps && mach <= mach_0 + math_utils::eps) ||
		(mach <= mach_1 + math_utils::eps && mach >= mach_0 - math_utils::eps)) {

		double delta = abs(last_state.getX(y + sign, x) - last_state.getX(y + sign, x - 1));
		
		if (abs(mach_1 - mach_0) < math_utils::eps) {
			state0 = pState1;

			x0 = last_state.getX(y + sign, x);
			y0 = last_state.getY(y - sign, x);
		}
		else {
			double h = (mach - mach_0) * delta / (mach_1 - mach_0);

			state0 = math_utils::interpolateState(pState1, pState0, delta, h);

			x0 = last_state.getX(y + sign, x) + h * sign;
			y0 = last_state.getY(y - sign, x);
		}
	}
	else if ((mach >= mach_1 - math_utils::eps && mach <= mach_2 + math_utils::eps) ||
		(mach <= mach_1 + math_utils::eps && mach >= mach_2 - math_utils::eps)) {

		if (abs(mach_1 - mach_2) < math_utils::eps) {
			state0 = pState1;

			x0 = last_state.getX(y + sign, x);
			y0 = last_state.getY(y - sign, x);
		}
		else {
			double delta = abs(last_state.getX(y + sign, x + 1) - last_state.getX(y + sign, x));
			double h = (mach - mach_1) * delta / (mach_2 - mach_1);

			state0 = math_utils::interpolateState(pState2, pState1, delta, h);

			x0 = last_state.getX(y + sign, x) + h * sign;
			y0 = last_state.getY(y - sign, x);
		}
	}
	else if ((mach <= mach_0 + math_utils::eps && mach_2 >= mach_0 - math_utils::eps) ||
		(mach >= mach_0 - math_utils::eps && mach_2 <= mach_0 + math_utils::eps)) {
		NodeState pState3 = last_state.getStateWithBoundary(y, x-1);
		double mach_3 = flow_utils::calcMach(pState3);

		if (abs(mach_3 - mach_0) < math_utils::eps) {
			state0 = pState1;

			x0 = last_state.getX(y + sign, x + 1);
			y0 = last_state.getY(y - sign, x);
		}
		else {
			double delta = abs(last_state.getY(y, x - 1) - last_state.getY(y + sign, x - 1));
			double h = (mach - mach_3) * delta / (mach_0 - mach_3);

			state0 = math_utils::interpolateState(pState3, pState0, delta, h);

			x0 = last_state.getX(y, x + 1);
			y0 = last_state.getY(y - sign, x + 1) + h * sign;
		}
	}
	else if ((mach >= mach_2 - math_utils::eps && mach_2 >= mach_0 - math_utils::eps) ||
		(mach <= mach_2 + math_utils::eps && mach_2 <= mach_0 + math_utils::eps)) {
		NodeState pState3 = last_state.getStateWithBoundary(y, x + 1);
		double mach_3 = flow_utils::calcMach(pState3);

		if (abs(mach_3 - mach_2) < math_utils::eps) {
			state0 = pState1;

			x0 = last_state.getX(y + sign, x - 1);
			y0 = last_state.getY(y - sign, x);
		}
		else {
			double delta = abs(last_state.getY(y, x + 1) - last_state.getY(y + sign, x + 1));
			double h = (mach - mach_3) * delta / (mach_2 - mach_3);

			state0 = math_utils::interpolateState(pState3, pState2, delta, h);

			x0 = last_state.getX(y + sign, x - 1);
			y0 = last_state.getY(y, x - 1) - h * sign;
		}
	}
	else {
		assert(0);
		std::cout << "Error: Can´t find line of constant Mach" << std::endl;
	}
}

void MachNeumannBoundary::calcExternalBoundary(State& last_state) {

	int y = last_state.getYSize() - 2;
	while (last_state.isWall(y, last_state.getXSize() / 2)) {
		y--;
	}
	int sx = 2;
	while (last_state.isWall(y, sx)) {
		sx++;
	}
	sx++;
	int ex = last_state.getXSize() - 2;
	while (last_state.isWall(y, ex)) {
		ex--;
	}
	ex--;
	double nxt_y = last_state.getY(y + 1, sx +1);
	NodeState& state = last_state[y][sx];
	auto boundary = state.boundary_;

	last_state[y + 1][sx] = last_state[y][sx];
	last_state[y + 1][sx].boundary_ = NULL;
	last_state[y + 1][ex] = last_state[y][ex];
	last_state[y + 1][ex].boundary_ = NULL;

	if (state.boundary_ != NULL &&
		state.boundaryType_ == BoundaryType::MACH_NEUMANN_BOUNDARY) {

		int last_node = 1;
		double last_node_x = last_state.getX(y, sx);
		NodeState last_node_state = last_state[y][sx];
		for (int x = sx +1; x < ex; x++) {
			double x0, y0;
			NodeState new_state;

			Boundary* boundary = NULL;
			for (auto* iter : *state.boundary_) {
				if (iter->getType() == BoundaryType::MACH_NEUMANN_BOUNDARY) {
					boundary = iter;
					break;
				}
			}

			MachNeumannBoundary* mach_neumann_boundary =
				static_cast<MachNeumannBoundary*>(boundary);
			mach_neumann_boundary->calcStateValue(x, y, last_state, x0, y0, new_state);

			double nxt_x = x0;

			// We want values in the mesh nodes. We assume here that the value is constant.
			if (abs(y0 - nxt_y) > 1.0e-8) {

				double m = (x0 - last_state.getX(y, x)) /
					(y0 - last_state.getY(y, x));

				nxt_x = m * (nxt_y - last_state.getY(y, x));
				nxt_x += last_state.getX(y, x);

			}

			int nxt_node = last_node;

			while (nxt_node + 1 < ex + 1 &&
				last_state.getX(y, nxt_node + 1) < nxt_x + 1.0e-8) {
				nxt_node++;

				double delta = nxt_x - last_node_x;
				double h = nxt_x - last_state.getX(y + 1, nxt_node);

				NodeState state = math_utils::interpolateState(last_node_state, new_state,
					delta, h);

				last_state[y + 1][nxt_node] = state;
				last_state[y + 1][nxt_node].boundary_ = state.boundary_;

			}

			last_node = nxt_node;
			last_node_x = x0;
			last_node_state = new_state;

		}

		double nxt_x = last_state.getX(y + 1, ex);
		NodeState new_state = last_state[y + 1][ex];
		for (int nxt_node = last_node + 1; nxt_node < ex; nxt_node++) {
			double delta = nxt_x - last_node_x;
			double h = nxt_x - last_state.getX(y + 1, nxt_node);

			NodeState state = math_utils::interpolateState(last_node_state, new_state,
				delta, h);

			last_state[y + 1][nxt_node] = state;
			last_state[y + 1][nxt_node].boundary_ = state.boundary_;
		}

	}

}