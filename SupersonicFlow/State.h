#pragma once
#include <vector>
#include "Mesh.h"
#include <memory>
#include "Boundary.h"
#include <math.h>
#include <algorithm>

class NodeState {
public:

	NodeState() {
		boundary_ = NULL;
	}
	NodeState(std::vector<double> v) {
		vals = v;
		boundary_ = NULL;
	}
	NodeState(int s) {
		vals.resize(s, 0.0);
		boundary_ = NULL;
	}

	void operator=(std::vector<double> v) {
		vals = v;
		boundary_ = NULL;
	}

	void operator=(NodeState state) {
		vals = state.vals;
		boundary_ = state.boundary_;
		done = state.done;
	}

	bool operator==(const NodeState& other) const {
		return other.vals == vals;
	}

	NodeState& operator*(double f) {
		for (int i = 0; i < vals.size(); i++) {
			vals[i] *= f;
		}
		return *this;
	}

	NodeState operator+(NodeState state) {
		NodeState st;
		for (int i = 0; i < vals.size(); i++) {
			st.add(vals[i] + state[i]);
		}
		return st;
	}

	NodeState operator-(NodeState state) const {
		NodeState st;
		for (int i = 0; i < vals.size(); i++) {
			st.add(vals[i] - state[i]);
		}
		return st;
	}

	NodeState operator/(double f) const {
		NodeState st;
		for (int i = 0; i < vals.size(); i++) {
			st.add(vals[i] / f);
		}
		return st;
	}

	void operator+=(NodeState state) {
		for (int i = 0; i < vals.size(); i++) {
			vals[i] += state[i];
		}
	}

	void operator-=(NodeState state) {
		for (int i = 0; i < vals.size(); i++) {
			vals[i] -= state[i];
		}
	}

	double& operator[](int i) {
		return vals[i];
	}

	int size() { return (int)vals.size(); }

	void add(double val) { vals.push_back(val); }

	void setBoundary(Boundary* boundary) {
		if (boundary_ == NULL) {
			boundary_ = std::shared_ptr<std::vector<Boundary*> >(new std::vector<Boundary*>());
		}
		boundary_->push_back(boundary);
		boundaryType_ = (BoundaryType)std::max((int)boundaryType_, (int)boundary->getType());
	}

	std::vector<double> vals;
	bool done = false;

	std::shared_ptr<std::vector<Boundary*> > boundary_;
	BoundaryType boundaryType_ = BoundaryType::EMPTY;
};

class State
{
public:
	State(Mesh & mesh);
	State();
	~State();

	static State createStateFrom(State& state);

	void initMesh(Mesh& mesh);
	Mesh* getMesh() { return mesh_; }

	std::vector<NodeState>& operator[](int i) { return state_[i]; }

	int getYSize() { return (int) state_.size(); }
	int getXSize() { return (int) state_[0].size(); }

	double getX(int y, int x) { return (*mesh_)[y][x].x; }
	double getY(int y, int x) { return (*mesh_)[y][x].y; }

	NodeState getStateWithBoundary(int y, int x);

	bool isWall(int y, int x) { return mesh_->isWall(y, x); }


	void exportTemperature(std::string file);
	void exportVelocityX(std::string file);
	void exportEnergy(std::string file);
	void exportPressure(std::string file);
	void exportVelocityY(std::string file);

private:
	std::vector<std::vector<NodeState> > state_;
	Mesh* mesh_;
};

NodeState operator*(double d, NodeState state);