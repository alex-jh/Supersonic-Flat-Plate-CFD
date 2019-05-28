#pragma once
#include "array2d_gpu.cuh"
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

__global__ void updatePredictedStepGpu(int imax, int jmax, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, double* delta_t, ArrayWrapper4D<double> U,
	ArrayWrapper4D<double> E, ArrayWrapper4D<double> F, ArrayWrapper4D<double> H, Array2DGpu<NODE_TYPE>* type, bool axilsymmetric,
	ArrayWrapper4D<double> U_predicted);

__global__ void updateCorrectedStepGpu(int imax, int jmax, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, double* delta_t, ArrayWrapper4D<double> U_predicted,
	ArrayWrapper4D<double> E, ArrayWrapper4D<double> F, ArrayWrapper4D<double> H, Array2DGpu<NODE_TYPE>* type, bool axilsymmetric,
	ArrayWrapper4D<double> U_corrected);

__global__ void encodeStateGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U);

__global__ void decodeStateGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U);

__global__ void calcEGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U,
	Array2DGpu<double>* tauxx, Array2DGpu<double>* tauxy, Array2DGpu<double>* qx, ArrayWrapper4D<double> E);

__global__ void calcFGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U,
	Array2DGpu<double>* tauyy, Array2DGpu<double>* tauxy, Array2DGpu<double>* qx, ArrayWrapper4D<double> F);

__global__ void calcHGpu(int imax, int jmax, Array2DGpu<double>* P, Array2DGpu<NODE_TYPE>* type,
	Array2DGpu<double>* tauh, ArrayWrapper4D<double> H);

__global__ void calcStresGpu(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* T, Array2DGpu<NODE_TYPE>* type,
	FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y, bool axilsymmetric,
	Array2DGpu<double>* tauxx, Array2DGpu<double>* tauyy, Array2DGpu<double>* tauxy);

__global__ void calcAxilsymmetricStressGpu(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* T, Array2DGpu<NODE_TYPE>* type,
	FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y, Array2DGpu<double>* tauh);

__global__ void calcHeatFluxGpu(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
	Array2DGpu<double>* T, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y,
	Array2DGpu<double>* qx, Array2DGpu<double>* qy);

class MacCormackSolverGpu
{
public:
	MacCormackSolverGpu(int imax, int jmax, bool axilsymmetric=false);
	~MacCormackSolverGpu();

	void Init(int imax, int jmax);
	void Destroy();

	//void Update(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
	//	Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
	//	Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type);

	void UpdatePredictor(double* delta_t, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, int imax, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type);

	void UpdateCorrector(double* delta_t, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, int imax, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type);

	void UpdatePredictedState(int imax, int jmax, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, double* delta_t, Array2DGpu<NODE_TYPE>* type);

	void UpdateCorrectedState(int imax, int jmax, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, double* delta_t, Array2DGpu<NODE_TYPE>* type);

	void EncodeState(int imax, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type);

	void DecodeState(int imax, int jmax, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double>& U);

	void CalcE(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	void CalcF(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	void CalcH(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
		Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
		Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	//void CalcStress(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
	//	Array2D<double>& u, Array2D<double>& v, Array2D<double>& T, Array2D<NODE_TYPE>& type,
	//	FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	//void CalcAxilsymmetricStress(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
	//	Array2D<double>& u, Array2D<double>& v, Array2D<double>& T, Array2D<NODE_TYPE>& type,
	//	FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

	//void CalcHeatFlux(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
	//	Array2D<double>& T, Array2D<NODE_TYPE>& type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y);

private:
	ArrayWrapper4D<double> U_;
	ArrayWrapper4D<double> U_predicted_;
	ArrayWrapper4D<double> U_corrected_;

	ArrayWrapper4D<double> E_;
	ArrayWrapper4D<double> F_;
	ArrayWrapper4D<double> H_;
	Array2DGpu<double>* tauxx_;
	Array2DGpu<double>* tauyy_;
	Array2DGpu<double>* tauxy_;
	Array2DGpu<double>* tauh_;
	Array2DGpu<double>* qx_;
	Array2DGpu<double>* qy_;

	bool axilsymmetric_;

public:
	int numThreads_;
	int numBlocks_;
};

