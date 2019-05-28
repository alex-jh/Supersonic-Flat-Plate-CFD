#pragma once
#include "array2d_gpu.cuh"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

template<class T>
bool WriteInFile(const Array2DGpu<T>* data, int imax, int jmax, const Array2DGpu<T>* deltax, const Array2DGpu<T>* deltay, std::string filename);