#pragma once
#include "math_utils.h"

class FlowParameters;

double CalcTStep(int imax, int jmax, double deltax, double deltay,
	const FlowParameters& params, const Array2D<double>& u, const Array2D<double>& v, 
	const Array2D<double>& rho, const Array2D<double>& P, const Array2D<double>& T, double K);

double CalcTStep2(int imax, int jmax, double deltax, double deltay,
	const FlowParameters& params, const Array2D<double>& u, const Array2D<double>& v,
	const Array2D<double>& rho, const Array2D<double>& P, const Array2D<double>& T, double K);
