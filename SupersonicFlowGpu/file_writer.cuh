#pragma once
#include "array2d_gpu.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

template<class T>
static bool WriteInFile(const Array2DGpu<T>* data, int imax, int jmax, double deltax, double deltay, std::string filename);