#include "pch.h"
#include "gtest/gtest.h"
#include "../SupersonicFlow/Boundary.h"
#include "../SupersonicFlow/MachNeumannBoundary.h"
#include "../SupersonicFlow/Mesh.h"
#include "../SupersonicFlow/State.h"
#include "../SupersonicFlow/Solver.h"

int size_x = 10;
int size_y = 10;

class BoundaryTest : public ::testing::Test {
public:
	BoundaryTest() {
		mesh_.createSimpleRectangularMesh(size_x, size_y);
		state_.initMesh(mesh_);
		state2_.initMesh(mesh_);
	}

protected:
	Mesh mesh_;
	State state_;

	State state2_;

friend class DirichletBoundaryTest;
friend class MachNeumannBoundaryTest;
};

class DirichletBoundaryTest : public BoundaryTest {
public:

	void SetUp() override {
	}

	void setLeftBoundary(NodeState& state_left_boundary) {
		for (int i = 0; i < state_.getYSize(); i++) {
			state_[i][0] = state_left_boundary;
		}

	}
};

TEST_F(DirichletBoundaryTest, ReturnsSameValue) {
	NodeState state_left_boundary;
	state_left_boundary.vals.push_back(1.0);
	setLeftBoundary(state_left_boundary);

	DirichletBoundary boundary;
	Solver solver;

	EXPECT_EQ(boundary.calcState(size_y/2, 0, state_, solver, state_), state_left_boundary);
}

class MachNeumannBoundaryTest : public BoundaryTest {
public:

	void SetUp() override {
		state_[1][1] = { 1, 1, 0, 1 };
		state_[1][2] = { 1, 1, 0, 5 };
		state_[1][3] = { 1, 1, 0, 17 };
		state_[2][1] = { 1, 1, 0, 1 };
		state_[2][2] = { 1, 1, 0, 10 };
		state_[2][3] = { 1, 1, 0, 26 };
		state_[2][4] = { 1, 1, 0, 37 };
		state_[2][1].boundary_ = &dirichlet_boundary;

		state_[1][7] = { 1, 1, 0, 1 };
		state_[1][8] = { 1, 1, 0, 5 };
		state_[1][9] = { 1, 1, 0, 17 };
		state_[2][6] = { 1, 1, 0, 1 };
		state_[2][7] = { 1, 1, 0, 6 };
		state_[2][8] = { 1, 1, 0, 10 };
		state_[2][9] = { 1, 1, 0, 26 };
		state_[2][9].boundary_ = &dirichlet_boundary;

		state_[9][9] = { 1, 1, 0, 1 };
		state_[9][8] = { 1, 1, 0, 5 };
		state_[9][7] = { 1, 1, 0, 17 };
		state_[8][9] = { 1, 1, 0, 1 };
		state_[8][8] = { 1, 1, 0, 6 };
		state_[8][7] = { 1, 1, 0, 10 };
		state_[8][6] = { 1, 1, 0, 26 };
		state_[8][9].boundary_ = &dirichlet_boundary;

		state_[9][1] = { 1, 1, 0, 1 };
		state_[9][2] = { 1, 1, 0, 17 };
		state_[9][3] = { 1, 1, 0, 24 };
		state_[8][1] = { 1, 1, 0, 1 };
		state_[8][2] = { 1, 1, 0, 6 };
		state_[8][3] = { 1, 1, 0, 10 };
		state_[8][4] = { 1, 1, 0, 26 };
		state_[8][1].boundary_ = &dirichlet_boundary;
		
		for (int i = 1; i <= 9; i++) {
			state2_[9][10 - i] = { 1, 1, 0, 2.0*i };
			state2_[8][10 - i] = { 1, 1, 0, 1.0*i };
		}
		state2_[9][1].boundary_ = &dirichlet_boundary;
		state2_[8][1].boundary_ = &dirichlet_boundary;
	}

	DirichletBoundary dirichlet_boundary;

};

TEST_F(MachNeumannBoundaryTest, BasicTestDirichletHorizontal) {

	MachNeumannBoundary boundary;
	Solver solver;

	std::pair<NodeState, double> pState = boundary.calcOutsideStateNodes(2, 1, state_);

	NodeState& state = pState.first;

	EXPECT_GT(state[3], state_[2][1][3]);
	EXPECT_LT(state[3], state_[2][3][3]);
}

TEST_F(MachNeumannBoundaryTest, BasicTestDirichletVertical) {

	MachNeumannBoundary boundary;
	Solver solver;

	std::pair<NodeState, double> pState = boundary.calcOutsideStateNodes(8, 1, state_);

	NodeState& state = pState.first;

	EXPECT_GT(state[3], state_[1][7][3]);
	EXPECT_LT(state[3], state_[2][7][3]);
	EXPECT_LT(pState.second, 1.0 - 1.0e-6);
	EXPECT_GT(pState.second, 1.0e-6);
}

TEST_F(MachNeumannBoundaryTest, BasicTestDirichletTop) {

	MachNeumannBoundary boundary;
	Solver solver;

	std::pair<NodeState, double> pState = boundary.calcOutsideStateNodes(8, 9, state_);

	NodeState& state = pState.first;

	EXPECT_GT(state[3], state_[8][7][3]);
	EXPECT_LT(state[3], state_[8][8][3]);
	EXPECT_LT(pState.second, 1.0 - 1.0e-6);
	EXPECT_GT(pState.second, 1.0e-6);
}

TEST_F(MachNeumannBoundaryTest, BasicTestDirichletVerticalTop) {

	MachNeumannBoundary boundary;
	Solver solver;

	std::pair<NodeState, double> pState = boundary.calcOutsideStateNodes(2, 9, state_);

	NodeState& state = pState.first;

	EXPECT_GT(state[3], state_[8][3][3]);
	EXPECT_GT(state[3], state_[9][3][3]);
	EXPECT_LT(pState.second, 1.0 - 1.0e-6);
	EXPECT_GT(pState.second, 1.0e-6);
}

TEST_F(MachNeumannBoundaryTest, Basic) {

	MachNeumannBoundary boundary;
	Solver solver;

	State new_state = State::createStateFrom(state2_);

	state2_[9][2].boundary_ = &boundary;


	boundary.calcExternalBoundary(state2_);

}