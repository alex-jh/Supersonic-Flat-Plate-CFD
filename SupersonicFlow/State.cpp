#include "State.h"
#include "Mesh.h"
#include "Boundary.h"
#include "Solver.h"
#include "FlowUtils.h"
#include <iostream>
#include <fstream>

NodeState operator*(double d, NodeState state) {
	return state * d;
}

State::State() {
}

State::State(Mesh& mesh) {
	initMesh(mesh);
}

State State::createStateFrom(State& state) {
	State new_state(*state.getMesh());
	return new_state;
}

void State::initMesh(Mesh& mesh) {
	mesh_ = &mesh;
	state_.resize(mesh.getYNodes());
	for (int i = 0; i < state_.size(); i++) {
		state_[i].resize(mesh.getXNodes());
	}
}


State::~State()
{
}

NodeState State::getStateWithBoundary(int y, int x) {

	if (x > 0 && x < getXSize() - 1 && y>0 && y < getYSize() - 1) {
		return state_[y][x];
	}

	int nxt_y = y;
	int nxt_x = x;

	if (x == 0) nxt_x = 1;
	else if (x == getXSize() - 1) nxt_x = getXSize() - 2;

	if (y == 0) nxt_y = 1;
	else if (y == getYSize() - 1) nxt_y = getYSize() - 2;

	switch (state_[nxt_y][nxt_x].boundaryType_)
	{
	case BoundaryType::DIRICHLET_BOUNDARY:
		return (*state_[nxt_y][nxt_x].boundary_)[0]->calcState(x, y, *this, Solver::GetSolver(), *this);
	case BoundaryType::UNKNOWN_BOUNDARY:
		return (*state_[nxt_y][nxt_x].boundary_)[0]->calcState(x, y, *this, Solver::GetSolver(), *this);
	default:
		std::cout << "Error: Unknown boundary type" << std::endl;
		return NodeState();
		break;
	}

}

void State::exportTemperature(std::string file) {

	std::ofstream outfile;
	outfile.open(file, std::ofstream::out);

	outfile << getYSize() << " " << getXSize() << std::endl;

	for (int i = 0; i < getXSize(); i++) {
		outfile << getX(0, i) << " ";
	}
	outfile << std::endl;
	for (int i = 0; i < getYSize(); i++) {
		outfile << getY(i, 0) << " ";
	}
	outfile << std::endl;

	for (int i = 0; i < getYSize(); i++) {
		for (int j = 0; j < getXSize(); j++) {
			double T = !isWall(i, j) ? flow_utils::calcTemperature(state_[i][j]) : 0;
			outfile << T << " ";
		}
		outfile << std::endl;
	}

	outfile.close();

}

void State::exportVelocityX(std::string file) {

	std::ofstream outfile;
	outfile.open(file, std::ofstream::out);

	outfile << getYSize() << " " << getXSize() << std::endl;

	for (int i = 0; i < getXSize(); i++) {
		outfile << getX(0, i) << " ";
	}
	outfile << std::endl;
	for (int i = 0; i < getYSize(); i++) {
		outfile << getY(i, 0) << " ";
	}
	outfile << std::endl;

	for (int i = 0; i < getYSize(); i++) {
		for (int j = 0; j < getXSize(); j++) {
			double T = !isWall(i, j) ? state_[i][j][1] / state_[i][j][0] : 0;
			outfile << T << " ";
		}
		outfile << std::endl;
	}

	outfile.close();

}

void State::exportEnergy(std::string file) {

	std::ofstream outfile;
	outfile.open(file, std::ofstream::out);

	outfile << getYSize() << " " << getXSize() << std::endl;

	for (int i = 0; i < getXSize(); i++) {
		outfile << getX(0, i) << " ";
	}
	outfile << std::endl;
	for (int i = 0; i < getYSize(); i++) {
		outfile << getY(i, 0) << " ";
	}
	outfile << std::endl;

	for (int i = 0; i < getYSize(); i++) {
		for (int j = 0; j < getXSize(); j++) {
			double T = !isWall(i, j) ? state_[i][j][3] / state_[i][j][0] : 0;
			outfile << T << " ";
		}
		outfile << std::endl;
	}

	outfile.close();

}

void State::exportVelocityY(std::string file) {

	std::ofstream outfile;
	outfile.open(file, std::ofstream::out);

	outfile << getYSize() << " " << getXSize() << std::endl;

	for (int i = 0; i < getXSize(); i++) {
		outfile << getX(0, i) << " ";
	}
	outfile << std::endl;
	for (int i = 0; i < getYSize(); i++) {
		outfile << getY(i, 0) << " ";
	}
	outfile << std::endl;

	for (int i = 0; i < getYSize(); i++) {
		for (int j = 0; j < getXSize(); j++) {
			double T = !isWall(i, j) ? state_[i][j][2] / state_[i][j][0] : 0;
			outfile << T << " ";
		}
		outfile << std::endl;
	}

	outfile.close();

}

void State::exportPressure(std::string file) {

	std::ofstream outfile;
	outfile.open(file, std::ofstream::out);

	outfile << getYSize() << " " << getXSize() << std::endl;

	for (int i = 0; i < getXSize(); i++) {
		outfile << getX(0, i) << " ";
	}
	outfile << std::endl;
	for (int i = 0; i < getYSize(); i++) {
		outfile << getY(i, 0) << " ";
	}
	outfile << std::endl;

	for (int i = 0; i < getYSize(); i++) {
		for (int j = 0; j < getXSize(); j++) {
			double T = !isWall(i, j) ? flow_utils::calcPressure(state_[i][j]) : 0;
			outfile << T << " ";
		}
		outfile << std::endl;
	}

	outfile.close();

}