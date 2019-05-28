#pragma once
#include "file_writer.h"
#include "maccormack_solver_gpu.h"

template
bool WriteInFile<double>(const Array2DGpu<double>* data, int imax, int jmax, const Array2DGpu<double>* deltax, const Array2DGpu<double>* deltay, std::string filename);

//template
//bool WriteInFile<NODE_TYPE>(const Array2DGpu<NODE_TYPE>* data, int imax, int jmax, const Array2DGpu<double>* deltax, const Array2DGpu<double>* deltay, std::string filename);

template<class T>
bool WriteInFile(const Array2DGpu<T>* data, int imax, int jmax, const Array2DGpu<T>* deltax, const Array2DGpu<T>* deltay, std::string filename) {

	T* buf;
	cudaMalloc((void**)&buf, imax * jmax * sizeof(T));

	copyArray2D<T> << <imax, jmax >> > (data, buf, imax*jmax);
	T* arr = new T[imax*jmax];
	cudaMemcpy(arr, buf, imax*jmax * sizeof(T), cudaMemcpyDeviceToHost);

	copyArray2D<T> << <2, 256 >> > (deltax, buf, 2 * imax);
	T* deltax_data = new T[2*imax];
	cudaMemcpy(deltax_data, buf, 2 * imax * sizeof(T), cudaMemcpyDeviceToHost);

	copyArray2D<T> << <2, 256 >> > (deltay, buf, 2 * jmax);
	T* deltay_data = new T[2*jmax];
	cudaMemcpy(deltay_data, buf, 2 * jmax * sizeof(T), cudaMemcpyDeviceToHost);

	std::ofstream file;
	file.open(filename + ".txt");
	file << std::fixed;
	file << imax << " " << jmax << std::endl;
	for (int i = 0; i < imax; i++) {
		for (int j = 0; j < jmax; j++) {
			file << std::setprecision(16) << *(deltax_data + imax + i) << " " <<
				*(deltay_data + jmax + j) << " " << *(arr + i * jmax + j) << std::endl;
		}
	}
	file.close();

	cudaFree(buf);

	delete deltay_data;
	delete deltax_data;
	delete arr;

	return true;
}