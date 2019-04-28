#include "maccormack.h"


MacCormack::MacCormack(int imax, int jmax) {
    for (int k = 0; k < 4; k++) {
        U_[k] = Array2D<double>(imax, jmax);
        U_predicted_[k] = Array2D<double>(imax, jmax);
        U_corrected_[k] = Array2D<double>(imax, jmax);
        E_[k] = Array2D<double>(imax, jmax);
        F_[k] = Array2D<double>(imax, jmax);
    }

    tauxx_ = Array2D<double>(imax, jmax);
    tauyy_ = Array2D<double>(imax, jmax);
    tauxy_ = Array2D<double>(imax, jmax);
    qx_ = Array2D<double>(imax, jmax);
    qy_ = Array2D<double>(imax, jmax);
}


MacCormack::~MacCormack()
{
}

void MacCormack::Update(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
    Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type) {

    UpdatePredictor(delta_t, delta_x, delta_y, imax, jmax, params, u, v, rho, P, T, e, type);

    UpdateCorrector(delta_t, delta_x, delta_y, imax, jmax, params, u, v, rho, P, T, e, type);
}

void MacCormack::UpdatePredictor(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
    Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type) {

    EncodeState(imax, jmax, params, u, v, rho, P, T, e);

    CalcE(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, REARWARDS_DIFFERENCES, CENTRAL_DIFFERENCES);
    CalcF(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, CENTRAL_DIFFERENCES, REARWARDS_DIFFERENCES);

    for (int i = 1; i < imax - 1; i++) {
        for (int j = 1; j < jmax - 1; j++) if (type.Get(i, j) == INSIDE) {
            for (int k = 0; k < 4; k++) {
                double x = E_[k].Get(i, j);
                double y = E_[k].Get(i + 1, j);
                double xx = F_[k].Get(i, j);
                double yy = F_[k].Get(i, j + 1);
                double a = E_[k].Get(i + 1, j) - E_[k].Get(i, j);
                double b = F_[k].Get(i, j + 1) - F_[k].Get(i, j);
                double tmp = -delta_t / delta_x * (E_[k].Get(i + 1, j) - E_[k].Get(i, j))
                    - delta_t / delta_y * (F_[k].Get(i, j + 1) - F_[k].Get(i, j));
                double tmp3 = U_predicted_[k].Get(i, j) + tmp;
                U_predicted_[k].Set(i, j, U_[k].Get(i, j) - delta_t / delta_x * (E_[k].Get(i + 1, j) - E_[k].Get(i, j))
                    - delta_t / delta_y * (F_[k].Get(i, j + 1) - F_[k].Get(i, j)));
            }
        }
    }

    DecodeState(imax, jmax, params, u, v, rho, P, T, e, U_predicted_);
}

void MacCormack::UpdateCorrector(double delta_t, double delta_x, double delta_y, int imax, int jmax, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
    Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type) {

    for (int i = 1; i < imax - 1; i++) {
        for (int j = 1; j < jmax - 1; j++) {
            for (int k = 0; k < 4; k++) {
                U_corrected_[k].Set(i, j, U_[k].Get(i, j));
            }
        }
    }

    EncodeState(imax, jmax, params, u, v, rho, P, T, e);

    CalcE(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, FORWARD_DIFFERENCES, CENTRAL_DIFFERENCES);
    CalcF(imax, jmax, delta_x, delta_y, params, u, v, rho, P, T, e, type, CENTRAL_DIFFERENCES, FORWARD_DIFFERENCES);

    for (int i = 1; i < imax - 1; i++) {
        for (int j = 1; j < jmax - 1; j++) if (type.Get(i, j) == INSIDE) {
            for (int k = 0; k < 4; k++) {
                double x = E_[k].Get(i, j);
                double y = E_[k].Get(i - 1, j);
                double xx = F_[k].Get(i, j);
                double yy = F_[k].Get(i, j - 1);
                double a = E_[k].Get(i, j) - E_[k].Get(i - 1, j);
                double b = F_[k].Get(i, j) - F_[k].Get(i, j - 1);
                double tmp = -delta_t / delta_x * (E_[k].Get(i, j) - E_[k].Get(i - 1, j))
                    - delta_t / delta_y * (F_[k].Get(i, j) - F_[k].Get(i, j - 1));
                double tmp2 = U_corrected_[k].Get(i, j);
                double tmp3 = U_predicted_[k].Get(i, j) + tmp;
                U_corrected_[k].Set(i, j, 0.5*(U_corrected_[k].Get(i, j) + U_predicted_[k].Get(i, j) - delta_t / delta_x * (E_[k].Get(i, j) - E_[k].Get(i - 1, j))
                    - delta_t / delta_y * (F_[k].Get(i, j) - F_[k].Get(i, j - 1))));
            }
        }
    }

    DecodeState(imax, jmax, params, u, v, rho, P, T, e, U_corrected_);
}

void MacCormack::EncodeState(int imax, int jmax, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
    Array2D<double>& T, Array2D<double>& e) {

    double rho_ij;
    double u_ij;
    double v_ij;

    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {

            rho_ij = rho.Get(i, j);
            u_ij = u.Get(i, j);
            v_ij = v.Get(i, j);

            U_[0].Set(i, j, rho_ij);
            U_[1].Set(i, j, rho_ij*u_ij);
            U_[2].Set(i, j, rho_ij*v_ij);
            U_[3].Set(i, j, rho_ij*(e.Get(i, j) + (u_ij*u_ij + v_ij * v_ij) / 2.));
        }
    }

}

void MacCormack::DecodeState(int imax, int jmax, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
    Array2D<double>& T, Array2D<double>& e, Array2D<double>(&U)[4]) {

    double rho_ij;
    double u_ij;
    double v_ij;
    double e_ij;

    for (int i = 1; i < imax - 1; i++) {
        for (int j = 1; j < jmax - 1; j++) {

            rho_ij = U[0].Get(i, j);
            u_ij = U[1].Get(i, j) / U[0].Get(i, j);
            v_ij = U[2].Get(i, j) / U[0].Get(i, j);
            e_ij = U[3].Get(i, j) / rho_ij - (u_ij*u_ij + v_ij * v_ij) / 2.0;

            u.Set(i, j, u_ij);
            v.Set(i, j, v_ij);
            rho.Set(i, j, rho_ij);
            T.Set(i, j, e_ij / params.cv);
            e.Set(i, j, e_ij);
            P.Set(i, j, params.R*rho_ij*T.Get(i, j));
        }
    }

}

void MacCormack::CalcE(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
    Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y) {

    double rho_ij;
    double u_ij;
    double v_ij;
    double p_ij;
    double tau_xx;
    double tau_xy;
    double q_x;

    CalcStress(imax, jmax, deltax, deltay, params, u, v, T, type, type_differences_x, type_differences_y);
    CalcHeatFlux(imax, jmax, deltax, deltay, params, T, type, type_differences_x, type_differences_y);

    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {

            rho_ij = rho.Get(i, j);
            u_ij = u.Get(i, j);
            v_ij = v.Get(i, j);
            p_ij = P.Get(i, j);
            tau_xx = tauxx_.Get(i, j);
            tau_xy = tauxy_.Get(i, j);
            q_x = qx_.Get(i, j);

            E_[0].Set(i, j, rho_ij*u_ij);
            E_[1].Set(i, j, rho_ij*u_ij*u_ij + p_ij - tau_xx);
            E_[2].Set(i, j, rho_ij*u_ij*v_ij - tau_xy);
            E_[3].Set(i, j, (U_[3].Get(i, j) + p_ij)*u_ij - u_ij * tau_xx - v_ij * tau_xy + q_x);
        }
    }

}

void MacCormack::CalcF(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
    Array2D<double>& T, Array2D<double>& e, Array2D<NODE_TYPE>& type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y) {

    double rho_ij;
    double u_ij;
    double v_ij;
    double p_ij;
    double tau_yy;
    double tau_xy;
    double q_y;

    CalcStress(imax, jmax, deltax, deltay, params, u, v, T, type, type_differences_x, type_differences_y);
    CalcHeatFlux(imax, jmax, deltax, deltay, params, T, type, type_differences_x, type_differences_y);

    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {

            rho_ij = rho.Get(i, j);
            u_ij = u.Get(i, j);
            v_ij = v.Get(i, j);
            p_ij = P.Get(i, j);
            tau_yy = tauyy_.Get(i, j);
            tau_xy = tauxy_.Get(i, j);
            q_y = qy_.Get(i, j);

            F_[0].Set(i, j, rho_ij*v_ij);
            F_[1].Set(i, j, rho_ij*u_ij*v_ij - tau_xy);
            F_[2].Set(i, j, rho_ij*v_ij*v_ij + p_ij - tau_yy);
            F_[3].Set(i, j, (U_[3].Get(i, j) + p_ij)*v_ij - u_ij * tau_xy - v_ij * tau_yy + q_y);
        }
    }

}

void MacCormack::CalcStress(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
    Array2D<double>& u, Array2D<double>& v, Array2D<double>& T, Array2D<NODE_TYPE>& type,
    FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y) {

    double u_ij;
    double v_ij;
    double t_ij;
    double ux_ij;
    double uy_ij;
    double vx_ij;
    double vy_ij;
    double mu_ij;

    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {

            u_ij = u.Get(i, j);
            v_ij = v.Get(i, j);
            t_ij = T.Get(i, j);
            mu_ij = ViscositySutherlandLaw(params, T.Get(i, j));

			if (i == 0 || (type.Get(i, j) == BOUNDARY && type.Get(i - 1, j) == OUTSIDE)) {
                ux_ij = (u.Get(i + 1, j) - u_ij) / deltax;
                vx_ij = (v.Get(i + 1, j) - v_ij) / deltax;
            }
            else if (i == imax - 1 || (type.Get(i, j) == BOUNDARY && type.Get(i + 1, j) == OUTSIDE)) {
                ux_ij = (u_ij - u.Get(i - 1, j)) / deltax;
                vx_ij = (v_ij - v.Get(i - 1, j)) / deltax;
            }
            else {
                switch (type_differences_x) {
                case FORWARD_DIFFERENCES:
                    ux_ij = (u.Get(i + 1, j) - u_ij) / deltax;
                    vx_ij = (v.Get(i + 1, j) - v_ij) / deltax;
                    break;
                case CENTRAL_DIFFERENCES:
                    ux_ij = (u.Get(i + 1, j) - u.Get(i - 1, j)) / (2 * deltax);
                    vx_ij = (v.Get(i + 1, j) - v.Get(i - 1, j)) / (2 * deltax);
                    break;
                case REARWARDS_DIFFERENCES:
                    ux_ij = (u_ij - u.Get(i - 1, j)) / deltax;
                    vx_ij = (v_ij - v.Get(i - 1, j)) / deltax;
                    break;
                }
            }

            if (j == 0 || (type.Get(i, j) == BOUNDARY && type.Get(i, j + 1) == OUTSIDE)) {
                uy_ij = (u.Get(i, j + 1) - u_ij) / deltay;
                vy_ij = (v.Get(i, j + 1) - v_ij) / deltay;
            }
            else if (j == jmax - 1 || (type.Get(i, j) == BOUNDARY && type.Get(i, j + 1) == OUTSIDE)) {
                uy_ij = (u_ij - u.Get(i, j - 1)) / deltay;
                vy_ij = (v_ij - v.Get(i, j - 1)) / deltay;
            }
            else {
                switch (type_differences_y) {
                case FORWARD_DIFFERENCES:
                    uy_ij = (u.Get(i, j + 1) - u_ij) / deltay;
                    vy_ij = (v.Get(i, j + 1) - v_ij) / deltay;
                    break;
                case CENTRAL_DIFFERENCES:
                    uy_ij = (u.Get(i, j + 1) - u.Get(i, j - 1)) / (2 * deltay);
                    vy_ij = (v.Get(i, j + 1) - v.Get(i, j - 1)) / (2 * deltay);
                    break;
                case REARWARDS_DIFFERENCES:
                    uy_ij = (u_ij - u.Get(i, j - 1)) / deltay;
                    vy_ij = (v_ij - v.Get(i, j - 1)) / deltay;
                    break;
                }
            }

            tauxx_.Set(i, j, -2. / 3 * mu_ij * (ux_ij + vy_ij) + 2 * mu_ij * ux_ij);
            tauyy_.Set(i, j, -2. / 3 * mu_ij * (ux_ij + vy_ij) + 2 * mu_ij * vy_ij);
            tauxy_.Set(i, j, mu_ij * (uy_ij + vx_ij));
        }
    }

}

void MacCormack::CalcHeatFlux(int imax, int jmax, double deltax, double deltay, const FlowParameters& params,
    Array2D<double>& T, Array2D<NODE_TYPE>& type, FiniteDifferencesType type_differences_x, FiniteDifferencesType type_differences_y) {

    double t_ij;
    double tx_ij;
    double ty_ij;
    double mu_ij;
    double k_ij;

    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {

            t_ij = T.Get(i, j);

            if (i == 0 || (type.Get(i, j) == BOUNDARY && type.Get(i - 1, j) == OUTSIDE)) {
                tx_ij = (T.Get(i + 1, j) - t_ij) / deltax;
            }
            else if (i == imax - 1 || (type.Get(i, j) == BOUNDARY && type.Get(i + 1, j) == OUTSIDE)) {
                tx_ij = (t_ij - T.Get(i - 1, j)) / deltax;
            }
            else {
                switch (type_differences_x) {
                case FORWARD_DIFFERENCES:
                    tx_ij = (T.Get(i + 1, j) - t_ij) / deltax;
                    break;
                case CENTRAL_DIFFERENCES:
                    tx_ij = (T.Get(i + 1, j) - T.Get(i - 1, j)) / (2 * deltax);
                    break;
                case REARWARDS_DIFFERENCES:
                    tx_ij = (t_ij - T.Get(i - 1, j)) / deltax;
                    break;
                }
            }

            if (j == 0 || (type.Get(i, j) == BOUNDARY && type.Get(i, j - 1) == OUTSIDE)) {
                ty_ij = (T.Get(i, j + 1) - t_ij) / deltay;
            }
            else if (j == jmax - 1 || (type.Get(i, j) == BOUNDARY && type.Get(i, j + 1) == OUTSIDE)) {
                ty_ij = (t_ij - T.Get(i, j - 1)) / deltay;
            }
            else {
                switch (type_differences_y) {
                case FORWARD_DIFFERENCES:
                    ty_ij = (T.Get(i, j + 1) - t_ij) / deltay;
                    break;
                case CENTRAL_DIFFERENCES:
                    ty_ij = (T.Get(i, j + 1) - T.Get(i, j - 1)) / (2 * deltay);
                    break;
                case REARWARDS_DIFFERENCES:
                    ty_ij = (t_ij - T.Get(i, j - 1)) / deltay;
                    break;
                }
            }

            mu_ij = ViscositySutherlandLaw(params, T.Get(i, j));
            k_ij = mu_ij * params.cp / params.Pr;

            qx_.Set(i, j, -k_ij * tx_ij);
            qy_.Set(i, j, -k_ij * ty_ij);
        }
    }

}