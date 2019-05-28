#include "maccormack_solver_gpu.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

void MacCormackSolverGpu::Init(int imax, int jmax) {

	for (int k = 0; k < 4; k++) {
		initArray2D<double> << <1, 1 >> > ((void*)U_[k], imax, jmax);
		//initArray2D<double> << <1, 1 >> > ((void*)(U2_.arr[k]), imax, jmax);
		initArray2D<double> << <1, 1 >> > ((void*)U_predicted_[k], imax, jmax);
		initArray2D<double> << <1, 1 >> > ((void*)U_corrected_[k], imax, jmax);

		initArray2D<double> << <1, 1 >> > ((void*)E_[k], imax, jmax);
		initArray2D<double> << <1, 1 >> > ((void*)F_[k], imax, jmax);
		initArray2D<double> << <1, 1 >> > ((void*)H_[k], imax, jmax);
		cudaDeviceSynchronize();
	}

	initArray2D<double> << <1, 1 >> > ((void*)tauxx_, imax, jmax);
	initArray2D<double> << <1, 1 >> > ((void*)tauyy_, imax, jmax);
	initArray2D<double> << <1, 1 >> > ((void*)tauxy_, imax, jmax);
	initArray2D<double> << <1, 1 >> > ((void*)tauh_, imax, jmax);
	initArray2D<double> << <1, 1 >> > ((void*)qx_, imax, jmax);
	initArray2D<double> << <1, 1 >> > ((void*)qy_, imax, jmax);

	cudaDeviceSynchronize();
}

void MacCormackSolverGpu::Destroy() {

	for (int k = 0; k < 4; k++) {
		destroyArray2D<double> << <1, 1 >> > (U_[k]);
		destroyArray2D<double> << <1, 1 >> > (U_predicted_[k]);
		destroyArray2D<double> << <1, 1 >> > (U_corrected_[k]);

		destroyArray2D<double> << <1, 1 >> > (E_[k]);
		destroyArray2D<double> << <1, 1 >> > (F_[k]);
		destroyArray2D<double> << <1, 1 >> > (H_[k]);
	}

	destroyArray2D<double> << <1, 1 >> > (tauxx_);
	destroyArray2D<double> << <1, 1 >> > (tauyy_);
	destroyArray2D<double> << <1, 1 >> > (tauxy_);
	destroyArray2D<double> << <1, 1 >> > (tauh_);
	destroyArray2D<double> << <1, 1 >> > (qx_);
	destroyArray2D<double> << <1, 1 >> > (qy_);

	cudaDeviceSynchronize();
}

__global__ void encodeStateGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {
			double rho_ij = rho->Get(i, j);
			double u_ij = u->Get(i, j);
			double v_ij = v->Get(i, j);

			U[0]->Set(i, j, rho_ij);
			U[1]->Set(i, j, rho_ij*u_ij);
			U[2]->Set(i, j, rho_ij*v_ij);
			U[3]->Set(i, j, rho_ij*(e->Get(i, j) + (u_ij*u_ij + v_ij * v_ij) / 2.));
		}

		idx += gridSize;
	}
}

void MacCormackSolverGpu::EncodeState(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type) {

	dim3 dimBlock(numThreads_, 1, 1);
	dim3 dimGrid(numBlocks_, 1, 1);

	encodeStateGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, params, u, v, rho, P, T, e, type, U_);

}

__global__ void decodeStateGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {
			double rho_ij = U[0]->Get(i, j);
			double u_ij = U[1]->Get(i, j) / U[0]->Get(i, j);
			double v_ij = U[2]->Get(i, j) / U[0]->Get(i, j);
			double e_ij = U[3]->Get(i, j) / rho_ij - (u_ij*u_ij + v_ij * v_ij) / 2.0;

			u->Set(i, j, u_ij);
			v->Set(i, j, v_ij);
			rho->Set(i, j, rho_ij);
			T->Set(i, j, e_ij / params->cv);
			e->Set(i, j, e_ij);
			P->Set(i, j, params->R*rho_ij*T->Get(i, j));
		}

		idx += gridSize;
	}
}

void MacCormackSolverGpu::DecodeState(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double>& U) {

	dim3 dimBlock(numThreads_, 1, 1);
	dim3 dimGrid(numBlocks_, 1, 1);

	decodeStateGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, params, u, v, rho, P, T, e, type, U);
}

__global__ void calcEGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U,
	Array2DGpu<double>* tauxx, Array2DGpu<double>* tauxy, Array2DGpu<double>* qx, ArrayWrapper4D<double> E) {
	
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {

			double rho_ij = rho->Get(i, j);
			double u_ij = u->Get(i, j);
			double v_ij = v->Get(i, j);
			double p_ij = P->Get(i, j);
			double tau_xx = tauxx->Get(i, j);
			double tau_xy = tauxy->Get(i, j);

			E[0]->Set(i, j, rho_ij*u_ij);
			E[1]->Set(i, j, rho_ij*u_ij*u_ij + p_ij - tau_xx);
			E[2]->Set(i, j, rho_ij*u_ij*v_ij - tau_xy);
			E[3]->Set(i, j, (U[3]->Get(i, j) + p_ij)*u_ij - u_ij * tau_xx - v_ij * tau_xy + qx->Get(i, j));
		}

		idx += gridSize;
	}
}

void MacCormackSolverGpu::CalcE(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y) {

	dim3 dimBlock(numThreads_, 1, 1);
	dim3 dimGrid(numBlocks_, 1, 1);

	calcStresGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, deltax, deltay, params, u, v, T, type, 
		type_differences_x, type_differences_y, axilsymmetric_, tauxx_, tauyy_, tauxy_);
	calcHeatFluxGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, deltax, deltay, params, T, type, type_differences_x, type_differences_y, qx_, qy_);

	calcEGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, params, u, v, rho, P, T, e, type, U_, tauxx_, tauxy_, qx_, E_);

}

__global__ void calcFGpu(int imax, int jmax, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, ArrayWrapper4D<double> U,
	Array2DGpu<double>* tauyy, Array2DGpu<double>* tauxy, Array2DGpu<double>* qy, ArrayWrapper4D<double> F) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {

			double rho_ij = rho->Get(i, j);
			double u_ij = u->Get(i, j);
			double v_ij = v->Get(i, j);
			double p_ij = P->Get(i, j);
			double tau_yy = tauyy->Get(i, j);
			double tau_xy = tauxy->Get(i, j);

			F[0]->Set(i, j, rho_ij*v_ij);
			F[1]->Set(i, j, rho_ij*u_ij*v_ij - tau_xy);
			F[2]->Set(i, j, rho_ij*v_ij*v_ij + p_ij - tau_yy);
			F[3]->Set(i, j, (U[3]->Get(i, j) + p_ij)*v_ij - u_ij * tau_xy - v_ij * tau_yy + qy->Get(i, j));
		}

		idx += gridSize;
	}
}

void MacCormackSolverGpu::CalcF(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y) {

	dim3 dimBlock(numThreads_, 1, 1);
	dim3 dimGrid(numBlocks_, 1, 1);

	calcStresGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, deltax, deltay, params, u, v, T, type,
		type_differences_x, type_differences_y, axilsymmetric_, tauxx_, tauyy_, tauxy_);
	calcHeatFluxGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, deltax, deltay, params, T, type, type_differences_x, type_differences_y, qx_, qy_);

	calcFGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, params, u, v, rho, P, T, e, type, U_, tauyy_, tauxy_, qy_, F_);
}

__global__ void calcHGpu(int imax, int jmax, Array2DGpu<double>* P, Array2DGpu<NODE_TYPE>* type,
	Array2DGpu<double>* tauh, ArrayWrapper4D<double> H) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {

			H[0]->Set(i, j, 0.0);
			H[1]->Set(i, j, 0.0);
			H[2]->Set(i, j, P->Get(i, j) - tauh->Get(i, j));
			H[3]->Set(i, j, 0.0);
		}

		idx += gridSize;
	}
}

void MacCormackSolverGpu::CalcH(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* rho, Array2DGpu<double>* P,
	Array2DGpu<double>* T, Array2DGpu<double>* e, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y) {

	if (!axilsymmetric_) {
		return;
	}

	dim3 dimBlock(numThreads_, 1, 1);
	dim3 dimGrid(numBlocks_, 1, 1);

	calcAxilsymmetricStressGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, deltax, deltay, params, u, v, T, type, 
		type_differences_x, type_differences_y, tauh_);

	calcHGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, P, type, tauh_, H_);

}


__global__ void calcStresGpu(int imax, int jmax, Array2DGpu<double>* deltax_v, Array2DGpu<double>* deltay_v, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* T, Array2DGpu<NODE_TYPE>* type, 
	FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y, bool axilsymmetric,
	Array2DGpu<double>* tauxx, Array2DGpu<double>* tauyy, Array2DGpu<double>* tauxy) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {

			double deltax = deltax_v->Get(i);
			double deltay = deltay_v->Get(j);

			double u_ij = u->Get(i, j);
			double v_ij = v->Get(i, j);
			double t_ij = T->Get(i, j);
			// This should not matter, because in d/dr we are multiplying F by r.
			// So this is to avoid nan.
			double r_ij = j != 0 ? deltay * j : deltay;
			double mu_ij = ViscositySutherlandLaw(*params, t_ij);
			double ux_ij, uy_ij, vx_ij, vy_ij;

			if (i == 0 || (type->Get(i, j) == BOUNDARY && type->Get(i - 1, j) == OUTSIDE)) {
				ux_ij = (u->Get(i + 1, j) - u_ij) / deltax;
				vx_ij = (v->Get(i + 1, j) - v_ij) / deltax;
			}
			else if (i == imax - 1 || (type->Get(i, j) == BOUNDARY && type->Get(i + 1, j) == OUTSIDE)) {
				ux_ij = (u_ij - u->Get(i - 1, j)) / deltax;
				vx_ij = (v_ij - v->Get(i - 1, j)) / deltax;
			}
			else {
				switch (type_differences_x) {
				case FORWARD_DIFFERENCES:
					ux_ij = (u->Get(i + 1, j) - u_ij) / deltax;
					vx_ij = (v->Get(i + 1, j) - v_ij) / deltax;
					break;
				case CENTRAL_DIFFERENCES:
					ux_ij = (u->Get(i + 1, j) - u->Get(i - 1, j)) / (2 * deltax);
					vx_ij = (v->Get(i + 1, j) - v->Get(i - 1, j)) / (2 * deltax);
					break;
				case REARWARDS_DIFFERENCES:
					ux_ij = (u_ij - u->Get(i - 1, j)) / deltax;
					vx_ij = (v_ij - v->Get(i - 1, j)) / deltax;
					break;
				}
			}

			if (j == 0 || (type->Get(i, j) == BOUNDARY && type->Get(i, j - 1) == OUTSIDE)) {
				uy_ij = (u->Get(i, j + 1) - u_ij) / deltay;
				vy_ij = (v->Get(i, j + 1) - v_ij) / deltay;
			}
			else if (j == jmax - 1 || (type->Get(i, j) == BOUNDARY && type->Get(i, j + 1) == OUTSIDE)) {
				uy_ij = (u_ij - u->Get(i, j - 1)) / deltay;
				vy_ij = (v_ij - v->Get(i, j - 1)) / deltay;
			}
			else {
				switch (type_differences_y) {
				case FORWARD_DIFFERENCES:
					uy_ij = (u->Get(i, j + 1) - u_ij) / deltay;
					vy_ij = (v->Get(i, j + 1) - v_ij) / deltay;
					break;
				case CENTRAL_DIFFERENCES:
					uy_ij = (u->Get(i, j + 1) - u->Get(i, j - 1)) / (2 * deltay);
					vy_ij = (v->Get(i, j + 1) - v->Get(i, j - 1)) / (2 * deltay);
					break;
				case REARWARDS_DIFFERENCES:
					uy_ij = (u_ij - u->Get(i, j - 1)) / deltay;
					vy_ij = (v_ij - v->Get(i, j - 1)) / deltay;
					break;
				}
			}

			if (!axilsymmetric) {
				tauxx->Set(i, j, -2. / 3 * mu_ij * (ux_ij + vy_ij) + 2 * mu_ij * ux_ij);
				tauyy->Set(i, j, -2. / 3 * mu_ij * (ux_ij + vy_ij) + 2 * mu_ij * vy_ij);
				tauxy->Set(i, j, mu_ij * (uy_ij + vx_ij));
			}
			else {
				tauxx->Set(i, j, -2. / 3 * mu_ij * (ux_ij + vy_ij + v_ij / r_ij) + 2 * mu_ij * ux_ij);
				tauyy->Set(i, j, -2. / 3 * mu_ij * (ux_ij + vy_ij + v_ij / r_ij) + 2 * mu_ij * vy_ij);
				tauxy->Set(i, j, mu_ij * (uy_ij + vx_ij));
			}
		}

		idx += gridSize;
	}

}

__global__ void calcAxilsymmetricStressGpu(int imax, int jmax, Array2DGpu<double>* deltax_v, Array2DGpu<double>* deltay_v, const FlowParameters* params,
	Array2DGpu<double>* u, Array2DGpu<double>* v, Array2DGpu<double>* T, Array2DGpu<NODE_TYPE>* type,
	FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y, Array2DGpu<double>* tauh) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {

			double deltax = deltax_v->Get(i);
			double deltay = deltay_v->Get(j);

			double u_ij = u->Get(i, j);
			double v_ij = v->Get(i, j);
			double t_ij = T->Get(i, j);
			double r_ij = deltay * j;
			double mu_ij = ViscositySutherlandLaw(*params, t_ij);
			double ux_ij, vy_ij;

			if (i == 0 || (type->Get(i, j) == BOUNDARY && type->Get(i - 1, j) == OUTSIDE)) {
				ux_ij = (u->Get(i + 1, j) - u_ij) / deltax;
			}
			else if (i == imax - 1 || (type->Get(i, j) == BOUNDARY && type->Get(i + 1, j) == OUTSIDE)) {
				ux_ij = (u_ij - u->Get(i - 1, j)) / deltax;
			}
			else {
				switch (type_differences_x) {
				case FORWARD_DIFFERENCES:
					ux_ij = (u->Get(i + 1, j) - u_ij) / deltax;
					break;
				case CENTRAL_DIFFERENCES:
					ux_ij = (u->Get(i + 1, j) - u->Get(i - 1, j)) / (2 * deltax);
					break;
				case REARWARDS_DIFFERENCES:
					ux_ij = (u_ij - u->Get(i - 1, j)) / deltax;
					break;
				}
			}

			if (j == 0 || (type->Get(i, j) == BOUNDARY && type->Get(i, j - 1) == OUTSIDE)) {
				vy_ij = (v->Get(i, j + 1) - v_ij) / deltay;
			}
			else if (j == jmax - 1 || (type->Get(i, j) == BOUNDARY && type->Get(i, j + 1) == OUTSIDE)) {
				vy_ij = (v_ij - v->Get(i, j - 1)) / deltay;
			}
			else {
				switch (type_differences_y) {
				case FORWARD_DIFFERENCES:
					vy_ij = (v->Get(i, j + 1) - v_ij) / deltay;
					break;
				case CENTRAL_DIFFERENCES:
					vy_ij = (v->Get(i, j + 1) - v->Get(i, j - 1)) / (2 * deltay);
					break;
				case REARWARDS_DIFFERENCES:
					vy_ij = (v_ij - v->Get(i, j - 1)) / deltay;
					break;
				}
			}

			tauh->Set(i, j, -2. / 3 * mu_ij * (ux_ij + vy_ij + v_ij / r_ij) + 2 * mu_ij * v_ij / r_ij);
		}

		idx += gridSize;
	}

}

__global__ void calcHeatFluxGpu(int imax, int jmax, Array2DGpu<double>* deltax_v, Array2DGpu<double>* deltay_v, const FlowParameters* params,
	Array2DGpu<double>* T, Array2DGpu<NODE_TYPE>* type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y,
	Array2DGpu<double>* qx, Array2DGpu<double>* qy) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < n) {
		unsigned int i = idx / jmax;
		unsigned int j = idx - i * jmax;

		if (type->Get(i, j) != OUTSIDE) {

			double deltax = deltax_v->Get(i);
			double deltay = deltay_v->Get(j);

			double t_ij = T->Get(i, j);
			double tx_ij, ty_ij;

			if (i == 0 || (type->Get(i, j) == BOUNDARY && type->Get(i - 1, j) == OUTSIDE)) {
				tx_ij = (T->Get(i + 1, j) - t_ij) / deltax;
			}
			else if (i == imax - 1 || (type->Get(i, j) == BOUNDARY && type->Get(i + 1, j) == OUTSIDE)) {
				tx_ij = (t_ij - T->Get(i - 1, j)) / deltax;
			}
			else {
				switch (type_differences_x) {
				case FORWARD_DIFFERENCES:
					tx_ij = (T->Get(i + 1, j) - t_ij) / deltax;
					break;
				case CENTRAL_DIFFERENCES:
					tx_ij = (T->Get(i + 1, j) - T->Get(i - 1, j)) / (2 * deltax);
					break;
				case REARWARDS_DIFFERENCES:
					tx_ij = (t_ij - T->Get(i - 1, j)) / deltax;
					break;
				}
			}

			if (j == 0 || (type->Get(i, j) == BOUNDARY && type->Get(i, j - 1) == OUTSIDE)) {
				ty_ij = (T->Get(i, j + 1) - t_ij) / deltay;
			}
			else if (j == jmax - 1 || (type->Get(i, j) == BOUNDARY && type->Get(i, j + 1) == OUTSIDE)) {
				ty_ij = (t_ij - T->Get(i, j - 1)) / deltay;
			}
			else {
				switch (type_differences_y) {
				case FORWARD_DIFFERENCES:
					ty_ij = (T->Get(i, j + 1) - t_ij) / deltay;
					break;
				case CENTRAL_DIFFERENCES:
					ty_ij = (T->Get(i, j + 1) - T->Get(i, j - 1)) / (2 * deltay);
					break;
				case REARWARDS_DIFFERENCES:
					ty_ij = (t_ij - T->Get(i, j - 1)) / deltay;
					break;
				}
			}

			double k_ij = ViscositySutherlandLaw(*params, t_ij) * params->cp / params->Pr;

			qx->Set(i, j, -k_ij * tx_ij);
			qy->Set(i, j, -k_ij * ty_ij);
		}
		idx += gridSize;
	}

}

__global__ void updatePredictedStepGpu(int imax, int jmax, Array2DGpu<double>* deltax_v, Array2DGpu<double>* deltay_v, double* delta_t, ArrayWrapper4D<double> U,
	ArrayWrapper4D<double> E, ArrayWrapper4D<double> F, ArrayWrapper4D<double> H, Array2DGpu<NODE_TYPE>* type, bool axilsymmetric,
	ArrayWrapper4D<double> U_predicted) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < 4 * n) {
		unsigned int k = idx / n;
		unsigned int i = (idx - k * n) / jmax;
		unsigned int j = (idx - k * n) - i * jmax;

		if (type->Get(i, j) == INSIDE) {
			double deltax = deltax_v->Get(i);
			double deltay = deltay_v->Get(j);

			double r = deltay_v->Get(jmax + j) * j;

			double delta = - (*delta_t) / deltax * (E[k]->Get(i + 1, j) - E[k]->Get(i, j))
				- (*delta_t) / deltay * (F[k]->Get(i, j + 1) - F[k]->Get(i, j));

			if (axilsymmetric) {
				delta = -(*delta_t) / deltax * (E[k]->Get(i + 1, j) - E[k]->Get(i, j))
					- (*delta_t) / deltay / r * ((r + deltay)*F[k]->Get(i, j + 1) - r * F[k]->Get(i, j)) +
					(*delta_t) * H[k]->Get(i, j) / r;
			}

			U_predicted[k]->Set(i, j, U[k]->Get(i, j) + delta);
		}

		idx += gridSize;
	}

}

void MacCormackSolverGpu::UpdatePredictedState(int imax, int jmax, Array2DGpu<double>* deltax, Array2DGpu<double>* deltay, double* delta_t, Array2DGpu<NODE_TYPE>* type) {

	dim3 dimBlock(numThreads_, 1, 1);
	dim3 dimGrid(numBlocks_, 1, 1);

	updatePredictedStepGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, deltax, deltay, delta_t, U_, E_, F_, H_, type, axilsymmetric_, U_predicted_);
}

__global__ void updateCorrectedStepGpu(int imax, int jmax, Array2DGpu<double>* deltax_v, Array2DGpu<double>* deltay_v, double* delta_t, ArrayWrapper4D<double> U_predicted,
	ArrayWrapper4D<double> E, ArrayWrapper4D<double> F, ArrayWrapper4D<double> H, Array2DGpu<NODE_TYPE>* type, bool axilsymmetric,
	ArrayWrapper4D<double> U_corrected) {

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	int n = imax * jmax;

	while (idx < 4 * n) {
		unsigned int k = idx / n;
		unsigned int i = (idx - k * n) / jmax;
		unsigned int j = (idx - k * n) - i * jmax;

		if (type->Get(i, j) == INSIDE) {

			double deltax = deltax_v->Get(i);
			double deltay = deltay_v->Get(j);

			double r = deltay_v->Get(jmax + j) * j;

			double delta = -(*delta_t) / deltax * (E[k]->Get(i, j) - E[k]->Get(i - 1, j))
				- (*delta_t) / deltay * (F[k]->Get(i, j) - F[k]->Get(i, j - 1));

			if (axilsymmetric) {
				delta = -(*delta_t) / deltax * (E[k]->Get(i, j) - E[k]->Get(i - 1, j))
					- (*delta_t) / deltay / r * (r*F[k]->Get(i, j) - (r - deltay)*F[k]->Get(i, j - 1)) +
					(*delta_t) * H[k]->Get(i, j) / r;
			}

			U_corrected[k]->Set(i, j, 0.5*(U_corrected[k]->Get(i, j) + U_predicted[k]->Get(i, j) + delta));
		}

		idx += gridSize;
	}

}

void MacCormackSolverGpu::UpdateCorrectedState(int imax, int jmax, Array2DGpu<double>* delta_x, Array2DGpu<double>* delta_y, double* delta_t, Array2DGpu<NODE_TYPE>* type) {

	dim3 dimBlock(numThreads_, 1, 1);
	dim3 dimGrid(numBlocks_, 1, 1);

	updateCorrectedStepGpu << <dimGrid, dimBlock, 1 >> > (imax, jmax, delta_x, delta_y, delta_t, U_predicted_, E_, F_, H_, type, axilsymmetric_, U_corrected_);
}