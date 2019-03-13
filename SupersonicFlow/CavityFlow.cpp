#include "CavityFlow.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include "Domain.h"
#include "Mesh.h"
#include "Boundary.h"
#include "MachNeumannBoundary.h"
#include "UnknownBoundary.h"
#include "WallBoundary.h"
#include "SymmetryBoundary.h"
#include "State.h"
#include "Solver.h"
#include "FlowUtils.h"
#include "CavityBoundaries.h"
#include <experimental/filesystem>
#include <stdio.h>
namespace fs = std::experimental::filesystem;

namespace cavity_flow {

	double u_c = 1;
	double P_c = 100e3;
	double r_c = 1.5;
	double m = 5000;
	double T_c = 302;
	double rho_0 = 0.01;

	void setInitialConditions(State& state) {

		flow_utils::FlowConstants::GetFlowConstants().T_c = T_c;
		flow_utils::FlowConstants::GetFlowConstants().P_c = P_c;
		flow_utils::FlowConstants::GetFlowConstants().u_c = u_c;
		flow_utils::FlowConstants::GetFlowConstants().r_c = r_c;
		flow_utils::FlowConstants::GetFlowConstants().Cp = 0.0;

		double T = 302 / T_c;
		double P = 100e3 / P_c;
		double u0 = 0 / u_c;
		double v0 = 0;

		double rho = flow_utils::calcDensity(P, T);
		double e0 = flow_utils::calcEnergy(P, T, u0, v0);

		NodeState initialState;
		initialState.add(rho);
		initialState.add(rho*u0);
		initialState.add(rho*v0);
		initialState.add(e0);

		for (int i = 0; i < state.getYSize(); i++) {
			for (int j = 0; j < state.getXSize(); j++) {
				if (!state.isWall(i, j)) {
					state[i][j] = initialState;
				}
			}
		}

	}

	NodeState getTopInitialCondition() {

		double T = 302 / T_c;
		double P = 100e3 / P_c;
		double u0 = u_c / u_c;
		double v0 = 0;

		double rho = flow_utils::calcDensity(P, T);
		double e0 = flow_utils::calcEnergy(P, T, u0, v0);

		NodeState initialState;
		initialState.add(rho);
		initialState.add(rho*u0);
		initialState.add(rho*v0);
		initialState.add(e0);

		return initialState;
	}

	CavityFlow::CavityFlow() {

		Domain domain("domain_square.bmp");

		mesh = Mesh("input_cavity.txt", domain);

		std::vector<std::vector<Node*> > left_boundaries = mesh.findLeftBoundary(domain, 2);
		std::vector<std::vector<Node*> > right_boundaries = mesh.findRightBoundary(domain, 2);
		std::vector<std::vector<Node*> > upper_boundaries = mesh.findUpperBoundary(domain, 2);
		std::vector<std::vector<Node*> > bottom_boundaries = mesh.findDownBoundary(domain, 2);

		boundaries.push_back(std::unique_ptr<Boundary>(new ConstantVelocityAndPressure()));
		ConstantVelocityAndPressure* boundary1 = static_cast<ConstantVelocityAndPressure*>(boundaries.back().get());

		boundaries.push_back(std::unique_ptr<Boundary>(new ConstantVelocityAndPressureX()));
		ConstantVelocityAndPressureX* boundary2 = static_cast<ConstantVelocityAndPressureX*>(boundaries.back().get());

		boundaries.push_back(std::unique_ptr<Boundary>(new ConstantVelocityAndPressureY()));
		ConstantVelocityAndPressureY* boundary3 = static_cast<ConstantVelocityAndPressureY*>(boundaries.back().get());

		state = State(mesh);

		setInitialConditions(state);

		{
			for (auto* it : upper_boundaries[0]) {
				// state[it->idxX][it->idxY].boundary_ = &mach_neumann_boundary;
				state[it->idxX][it->idxY].setBoundary(boundary1);
			}

			for (auto* it : bottom_boundaries[0]) {
				state[it->idxX][it->idxY].setBoundary(boundary3);
			}

			for (auto* it : left_boundaries[0]) {
				state[it->idxX][it->idxY].setBoundary(boundary2);
			}

			for (auto* it : right_boundaries[0]) {
				state[it->idxX][it->idxY].setBoundary(boundary2);
			}
		}

		NodeState initial_condition = getTopInitialCondition();

		{
			for (auto* it : upper_boundaries[0]) {
				auto boundary = state[it->idxX][it->idxY].boundary_;
				state[it->idxX][it->idxY] = initial_condition;
				state[it->idxX][it->idxY].boundary_ = boundary;
			}
		}
	}


	CavityFlow::~CavityFlow()
	{
	}

	void CavityFlow::Run() {
		Solver& solver = Solver::GetSolver();
		solver.axilsymmetric = false;
		solver.delta_t = .00001;

		int iter_log = 10;
		std::string path = "C:\\Users\\Alex\\Documents\\cavity_out\\";

		fs::create_directories(path);

		for (int i = 0; i <= 1000; i++) {
			std::cout << "Iteration: " << i << std::endl;
			state = solver.solve(state);

			printf("%d %.12f %.12f\n", i, state[97][3][3], flow_utils::calcPressure(state[97][3]));

			if (i % iter_log == 0) {
				std::ostringstream stringStream;
				stringStream << path << "temperature_" << i << ".txt";

				state.exportTemperature(stringStream.str());
			}

			if (i % iter_log == 0) {
				std::ostringstream stringStream;
				stringStream << path << "velocityX_" << i << ".txt";

				state.exportVelocityX(stringStream.str());
			}

			if (i % iter_log == 0) {
				std::ostringstream stringStream;
				stringStream << path << "velocityY_" << i << ".txt";

				state.exportVelocityY(stringStream.str());
			}

			if (i % iter_log == 0) {
				std::ostringstream stringStream;
				stringStream << path << "pressure_" << i << ".txt";

				state.exportPressure(stringStream.str());
			}

			if (i % iter_log == 0) {
				std::ostringstream stringStream;
				stringStream << path << "energy_" << i << ".txt";

				state.exportEnergy(stringStream.str());
			}
		}
	}
}