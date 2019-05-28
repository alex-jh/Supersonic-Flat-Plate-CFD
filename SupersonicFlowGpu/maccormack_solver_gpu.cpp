#include "maccormack_solver_gpu.h"
#include "file_writer.h"

MacCormackSolverGpu::MacCormackSolverGpu(int imax, int jmax, bool axylsymmetric) {
	for (int k = 0; k < 4; k++) {
		cudaMalloc((void**)&(U_.arr[k]), sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&(U_predicted_.arr[k]), sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&(U_corrected_.arr[k]), sizeof(Array2DGpu<double>));
		
		cudaMalloc((void**)&(E_.arr[k]), sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&(F_.arr[k]), sizeof(Array2DGpu<double>));
		cudaMalloc((void**)&(H_.arr[k]), sizeof(Array2DGpu<double>));
	}

	cudaMalloc((void**)&tauxx_, sizeof(Array2DGpu<double>));
	cudaMalloc((void**)&tauyy_, sizeof(Array2DGpu<double>));
	cudaMalloc((void**)&tauxy_, sizeof(Array2DGpu<double>));
	cudaMalloc((void**)&tauh_, sizeof(Array2DGpu<double>));
	cudaMalloc((void**)&qx_, sizeof(Array2DGpu<double>));
	cudaMalloc((void**)&qy_, sizeof(Array2DGpu<double>));

	axilsymmetric_ = axylsymmetric;

	Init(imax, jmax);
}


MacCormackSolverGpu::~MacCormackSolverGpu()
{
	Destroy();

	for (int k = 0; k < 4; k++) {
		cudaFree(U_[k]);
		cudaFree(U_predicted_[k]);
		cudaFree(U_corrected_[k]);
		cudaFree(E_[k]);
		cudaFree(F_[k]);
		cudaFree(H_[k]);
	}

	cudaFree(tauxx_);
	cudaFree(tauyy_);
	cudaFree(tauxy_);
	cudaFree(tauh_);
	cudaFree(qx_);
	cudaFree(qy_);
}

//void MacCormack::Update(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
//	Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
//	Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type) {
//
//	UpdatePredictor(delta_t, delta_x, delta_y, imax, jmax, params, u, v, rho, P, T, e, type);
//
//	UpdateCorrector(delta_t, delta_x, delta_y, imax, jmax, params, u, v, rho, P, T, e, type);
//}

void MacCormackSolverGpu::UpdatePredictor(double* delta_t, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type) {

	EncodeState(imax, jmax, params, u, v, rho, P, T, e, type);
	// cudaDeviceSynchronize();

	CalcE(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, REARWARDS_DIFFERENCES, CENTRAL_DIFFERENCES);
	CalcF(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, CENTRAL_DIFFERENCES, REARWARDS_DIFFERENCES);
	CalcH(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, CENTRAL_DIFFERENCES, REARWARDS_DIFFERENCES);
	// cudaDeviceSynchronize();

	UpdatePredictedState(imax, jmax, delta_x, delta_y, delta_t, type);
	// cudaDeviceSynchronize();

	DecodeState(imax, jmax, params, u, v, rho, P, T, e, type, U_predicted_);
	// cudaDeviceSynchronize();
}

void MacCormackSolverGpu::UpdateCorrector(double* delta_t, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type) {

	for (int k = 0; k < 4; k++) {
		U_corrected_[k]->copyArray2D(imax*jmax, U_[k], numThreads_, numBlocks_);
	}

	// cudaDeviceSynchronize();
	EncodeState(imax, jmax, params, u, v, rho, P, T, e, type);
	// cudaDeviceSynchronize();

	CalcE(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, FORWARD_DIFFERENCES, CENTRAL_DIFFERENCES);
	CalcF(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, CENTRAL_DIFFERENCES, FORWARD_DIFFERENCES);
	CalcH(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, CENTRAL_DIFFERENCES, FORWARD_DIFFERENCES);
	// cudaDeviceSynchronize();

	UpdateCorrectedState(imax, jmax, delta_x, delta_y, delta_t, type);
	// cudaDeviceSynchronize();

	DecodeState(imax, jmax, params, u, v, rho, P, T, e, type, U_corrected_);
	// cudaDeviceSynchronize();
}

