#include "flow_utils.h"
#include "flow_parameters.h"
#include <math.h>
#include <algorithm>

double ViscositySutherlandLaw(const FlowParameters& params, double T) {
    return params.mu*T / params.T_0*sqrt(T / params.T_0)*
		(params.T_0 + 110) / (T + 110);
}
