#include "SupersonicNozzle.h"
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
#include <experimental/filesystem>
#include <stdio.h>
namespace fs = std::experimental::filesystem;

namespace supersonic_nozzle {
	double u_j = 900;
	double P_j = 15e3;
	double r_j = 1.5;
	double T_j = 650;
	double u_c = u_j;
	double P_c = P_j;
	double r_c = 1.5;
	double T_c = T_j;
	double m = 5000;

	void setInitialConditions(State& state) {

		flow_utils::FlowConstants::GetFlowConstants().T_c = T_c;
		flow_utils::FlowConstants::GetFlowConstants().P_c = P_c;
		flow_utils::FlowConstants::GetFlowConstants().u_c = u_c;
		flow_utils::FlowConstants::GetFlowConstants().r_c = r_c;
		flow_utils::FlowConstants::GetFlowConstants().Cp = 1.005;

		double T = 302 / T_c;
		double P = 30e3 / P_c;
		double u0 = 0 / u_c;
		double v0 = 0;

		double rho = flow_utils::calcDensity(P, T);
		double e0 = flow_utils::calcEnergy(P, T, u0, v0);

		NodeState initialState;
		initialState.add(rho);
		initialState.add(rho*u0);
		initialState.add(rho*v0);
		initialState.add(rho*e0);

		for (int i = 0; i < state.getYSize(); i++) {
			for (int j = 0; j < state.getXSize(); j++) {
				if (!state.isWall(i, j)) {
					state[i][j] = initialState;
				}
			}
		}

	}

	NodeState getJetInitialCondition() {
		NodeState state;

		double u0 = u_j / u_c;
		double v0 = 0;
		double P = P_j / P_c;
		double A = r_c * r_c*acos(-1.0);
		double m = 5000;

		double T = T_j / T_c;

		double rho = flow_utils::calcDensity(P, T);
		double e0 = flow_utils::calcEnergy(P, T, u0, v0);

		state.add(rho);
		state.add(rho*u0);
		state.add(rho*v0);
		state.add(rho*e0);

		return state;
	}

	SupersonicNozzle::SupersonicNozzle() {
		Domain domain("domain.bmp");

		mesh = Mesh("input.txt", domain);

		std::vector<std::vector<Node*> > left_boundaries = mesh.findLeftBoundary(domain, 2);
		std::vector<std::vector<Node*> > right_boundaries = mesh.findRightBoundary(domain, 2);
		std::vector<std::vector<Node*> > upper_boundaries = mesh.findUpperBoundary(domain, 2);
		std::vector<std::vector<Node*> > bottom_boundaries = mesh.findDownBoundary(domain, 2);

		std::sort(bottom_boundaries.begin(), bottom_boundaries.end(),
			[](const std::vector<Node*> & a, const std::vector<Node*> & b) -> bool
		{
			return a.size() < b.size();
		});

		std::sort(upper_boundaries.begin(), upper_boundaries.end(),
			[](const std::vector<Node*> & a, const std::vector<Node*> & b) -> bool
		{
			return a.size() < b.size();
		});

		boundaries.push_back(std::unique_ptr<Boundary>(new DirichletBoundary()));
		DirichletBoundary* dirichlet_boundary = static_cast<DirichletBoundary*>(boundaries.back().get());

		boundaries.push_back(std::unique_ptr<Boundary>(new MachNeumannBoundary()));
		MachNeumannBoundary* mach_neumann_boundar = static_cast<MachNeumannBoundary*>(boundaries.back().get());

		boundaries.push_back(std::unique_ptr<Boundary>(new UnknownBoundary()));
		UnknownBoundary* unknown_boundary = static_cast<UnknownBoundary*>(boundaries.back().get());

		boundaries.push_back(std::unique_ptr<Boundary>(new WallBoundary()));
		WallBoundary* wall_boundary = static_cast<WallBoundary*>(boundaries.back().get());

		boundaries.push_back(std::unique_ptr<Boundary>(new SymmetryBoundary()));
		SymmetryBoundary* symmetry_boundary = static_cast<SymmetryBoundary*>(boundaries.back().get());

		state = State(mesh);

		setInitialConditions(state);

		{
			for (auto* it : upper_boundaries[1]) {
				// state[it->idxX][it->idxY].boundary_ = &mach_neumann_boundary;
				state[it->idxX][it->idxY].setBoundary(unknown_boundary);
			}

			for (auto* it : bottom_boundaries[1]) {
				state[it->idxX][it->idxY].setBoundary(symmetry_boundary);
			}
		}

		{
			for (auto* it : upper_boundaries[0]) {
				state[it->idxX][it->idxY].setBoundary(wall_boundary);
			}

			for (auto* it : bottom_boundaries[0]) {
				state[it->idxX][it->idxY].setBoundary(wall_boundary);
			}
		}

		{
			for (auto* it : left_boundaries[2]) {
				state[it->idxX][it->idxY].setBoundary(wall_boundary);
			}

			for (auto* it : left_boundaries[1]) {
				state[it->idxX][it->idxY].setBoundary(dirichlet_boundary);
			}
		}

		for (auto* it : right_boundaries[0]) {
			state[it->idxX][it->idxY].setBoundary(unknown_boundary);
		}

		NodeState jet_initial_condition = getJetInitialCondition();

		{
			for (auto* it : left_boundaries[0]) {
				std::shared_ptr<std::vector<Boundary*> > tmp = state[it->idxX][it->idxY].boundary_;
				state[it->idxX][it->idxY] = jet_initial_condition;
				state[it->idxX][it->idxY].boundary_ = tmp;
				state[it->idxX][it->idxY].setBoundary(dirichlet_boundary);
			}
		}
	}


	SupersonicNozzle::~SupersonicNozzle()
	{
	}

	void SupersonicNozzle::Run() {
		Solver& solver = Solver::GetSolver();
		solver.axilsymmetric = true;
		solver.delta_t = .0001;

		int iter_log = 1;
		std::string path = "C:\\Users\\Alex\\Documents\\cavity_out\\";

		fs::create_directories(path);

		for (int i = 0; i <= 10000; i++) {
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