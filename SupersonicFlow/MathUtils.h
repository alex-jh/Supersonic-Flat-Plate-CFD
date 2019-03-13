#pragma once
#include "State.h"
#include <assert.h>

namespace math_utils {

	const double eps = 1.0e-8;

	double interpolate(double a, double b, double delta, double h) {
		return (a*h + b * (delta - h)) / delta;
	}

	NodeState interpolateState(NodeState& st1, NodeState& st2, double delta, double h) {
		NodeState st;
		for (int i = 0; i < st1.size(); i++) {
			st.add(interpolate(st1[i], st2[i], delta, h));
		}

		return st;
	}

	NodeState parabolicInterpolateState(std::vector<NodeState>& vStates, std::vector<double>& vx, double x) {
		double l0 = (x - vx[1])*(x - vx[2]) /
			((vx[0] - vx[1])*(vx[0] - vx[2]));
		double l1 = (x - vx[0])*(x - vx[2]) /
			((vx[1] - vx[0])*(vx[1] - vx[2]));
		double l2 = (x - vx[0])*(x - vx[1]) /
			((vx[2] - vx[0])*(vx[2] - vx[1]));

		NodeState st;
		bool check = true;
		for (int i = 0; i < vStates[0].size(); i++) {
			double val = vStates[0][i] * l0 + vStates[1][i] * l1 +
				vStates[2][i] * l2;
			st.add(val);
			if ((val < vStates[0][i] && val < vStates[2][i]) ||
				(val > vStates[0][i] && val > vStates[2][i])) {
				check = false;
			}
		}

		if (!check) {
			if ((x >= vx[0] && x <= vx[1]) ||
				(x <= vx[0] && x >= vx[1])) {
				double delta = abs(vx[0] - vx[1]);
				double h = abs(x - vx[1]);
				return  interpolateState(vStates[0], vStates[1], delta, h);
			}
			else if ((x >= vx[2] && x <= vx[1]) ||
				(x <= vx[2] && x >= vx[1])) {
				double delta = abs(vx[2] - vx[1]);
				double h = abs(x - vx[2]);
				return  interpolateState(vStates[1], vStates[2], delta, h);
			}
			else {
				assert(0&&"Interpolation error");
				std::cout << "Interpolation error" << std::endl;
			}
		}

		return st;

	}

	double parabolicInterpolate(std::vector<double>& vStates, std::vector<double>& vx, double x) {
		double l0 = (x - vx[1])*(x - vx[2]) /
			((vx[0] - vx[1])*(vx[0] - vx[2]));
		double l1 = (x - vx[0])*(x - vx[2]) /
			((vx[1] - vx[0])*(vx[1] - vx[2]));
		double l2 = (x - vx[0])*(x - vx[1]) /
			((vx[2] - vx[0])*(vx[2] - vx[1]));

		double val = vStates[0] * l0 + vStates[1] * l1 +
			vStates[2] * l2;

		bool check = true;
		if ((val < vStates[0] && val < vStates[2]) ||
			(val > vStates[0] && val > vStates[2])) {
			check = false;
		}

		if (!check) {
			if ((x >= vx[0] && x <= vx[1]) ||
				(x <= vx[0] && x >= vx[1])) {
				double delta = abs(vx[0] - vx[1]);
				double h = abs(x - vx[1]);
				return  interpolate(vStates[0], vStates[1], delta, h);
			}
			else if ((x >= vx[2] && x <= vx[1]) ||
				(x <= vx[2] && x >= vx[1])) {
				double delta = abs(vx[2] - vx[1]);
				double h = abs(x - vx[2]);
				return  interpolate(vStates[1], vStates[2], delta, h);
			}
			else {
				assert(0 && "Interpolation error");
				std::cout << "Interpolation error" << std::endl;
			}
		}

		return val;

	}

}