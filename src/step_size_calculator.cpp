#include "step_size_calculator.h"
#include "flow_parameters.h"
#include <math.h>
#include "flow_utils.h"
#include <algorithm>

double CalcTStep(int imax, int jmax, double deltax, double deltay,
    const FlowParameters& params, const Array2D<double>& u, const Array2D<double>& v,
    const Array2D<double>& rho, const Array2D<double>& P, const Array2D<double>& T, double K) {

    double delta_t = 1.0;

    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {

            double u_ij = u.Get(i, j);
            double v_ij = v.Get(i, j);

            double delta_t_ij = abs(u_ij) / deltax;
            delta_t_ij += abs(v_ij) / deltay;

            double a_ij = sqrt(std::max(0.0, P.Get(i, j)*params.gamma / rho.Get(i, j)));
            delta_t_ij += a_ij * sqrt(1. / deltax / deltax + 1. / deltay / deltay);

            double mu_ij = ViscositySutherlandLaw(params, T.Get(i, j));
            double vv_ij = abs(4. / 3.*mu_ij*(params.gamma*mu_ij / params.Pr) / rho.Get(i, j));
            delta_t_ij += 2 * vv_ij*(1. / deltax / deltax + 1. / deltay / deltay);

            delta_t_ij = 1.0 / delta_t_ij;

            delta_t = std::min(delta_t, delta_t_ij*K);
        }
    }

    return delta_t;
}
