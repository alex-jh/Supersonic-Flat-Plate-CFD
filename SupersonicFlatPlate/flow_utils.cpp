#include "flow_utils.h"
#include "flow_parameters.h"
#include <math.h>
#include <algorithm>

double ViscositySutherlandLaw(const FlowParameters& params, double T) {
    return params.mu*pow(std::max(0.0, T / params.T_0), 1.5)*(params.T_0 + 110) / (T + 110);
}
