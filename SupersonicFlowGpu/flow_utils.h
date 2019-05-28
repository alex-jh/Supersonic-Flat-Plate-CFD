#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include "flow_parameters.h"
class FlowParameters;


//__device__ double ViscositySutherlandLaw(const FlowParameters& params, double T);

inline __device__ __host__ double ViscositySutherlandLaw(const FlowParameters& params, double T) {
	return params.mu*T / params.T_0*sqrt(T / params.T_0)*
		(params.T_0 + 110) / (T + 110);
}