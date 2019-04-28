#pragma once
#include "math_utils.h"
#include "flow_utils.h"
#include "flow_parameters.h"

enum FiniteDifferencesType {
	FORWARD_DIFFERENCES,
	CENTRAL_DIFFERENCES,
	REARWARDS_DIFFERENCES
};

typedef enum NodeType {
	INSIDE = 0,
	BOUNDARY,
	OUTSIDE
} NODE_TYPE;

class MacCormack
{
public:
	MacCormack(int imax, int jmax);
	~MacCormack();

	void Update(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type);

	void UpdatePredictor(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type);

	void UpdateCorrector(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type);

	void EncodeState(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e);

	void DecodeState(int imax, int jmax, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e, Array2D<double>(&U)[4]);

	void CalcE(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	void CalcF(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
		Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	void CalcStress(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
		Array2D<double>& u, Array2D<double>& v, Array2D<double>& T, Array2D<NODE_TYPE>& type,
		FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	void CalcHeatFlux(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
		Array2D<double>& T, Array2D<NODE_TYPE>& type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);
private:
	Array2D<double> U_[4];
	Array2D<double> U_predicted_[4];
	Array2D<double> U_corrected_[4];

	Array2D<double> E_[4];
	Array2D<double> F_[4];
	Array2D<double> tauxx_;
	Array2D<double> tauyy_;
	Array2D<double> tauxy_;
	Array2D<double> qx_;
	Array2D<double> qy_;
};

