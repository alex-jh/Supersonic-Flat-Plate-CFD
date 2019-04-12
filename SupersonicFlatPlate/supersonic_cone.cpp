#include "supersonic_cone.h"
#include "math.h"
#include "step_size_calculator.h"
#include <algorithm>
#include <iostream>
#include "file_writer.h"

namespace supersonic_cone {

    const int IMAX = 120;
    const int JMAX = 70;
    const double pi = acos(-1.0);
    const double angle = 0;
    const int margin = 0;
    const int height = 0;

    SupersonicCone::SupersonicCone() : imax_(IMAX), jmax_(JMAX), maccormack_solver_(IMAX, JMAX)
    {
        flow_parameters_.mu = 1.7894e-5;
        flow_parameters_.T_0 = 288.16;

        double M_inf = 3;
        double a_inf = 340.28;
        double u_inf = a_inf * M_inf;
        double plate_length = 0.00001;
        double R = 287;
        double P_inf = 101325;
        double T_inf = 288.16;
        double rho_inf = P_inf / T_inf / R;
        double mu_inf = ViscositySutherlandLaw(flow_parameters_, T_inf);
        double L = sqrt(mu_inf*plate_length / rho_inf / u_inf);

        flow_parameters_.M_inf = M_inf;
        flow_parameters_.plate_length = 0.00001 / L;
        flow_parameters_.a_inf = a_inf / u_inf;
        flow_parameters_.P_inf = P_inf / P_inf;
        flow_parameters_.T_inf = T_inf / T_inf;
        flow_parameters_.T_wall = flow_parameters_.T_inf;
        flow_parameters_.gamma = 1.4;
        flow_parameters_.R = 287 / u_inf / u_inf * T_inf;
        flow_parameters_.Pr = 0.71;
        flow_parameters_.cv = 0.718;
        flow_parameters_.cp = 1.01;
        flow_parameters_.mu /= (P_inf * (L / u_inf));
        flow_parameters_.T_0 /= T_inf;

        maxit = 20000;

        T_ = Array2D<double>(imax_, jmax_);
        u_ = Array2D<double>(imax_, jmax_);
        v_ = Array2D<double>(imax_, jmax_);
        P_ = Array2D<double>(imax_, jmax_);
        rho_ = Array2D<double>(imax_, jmax_);
        M_ = Array2D<double>(imax_, jmax_);
        e_ = Array2D<double>(imax_, jmax_);
        outside_ = Array2D<int>(imax_, jmax_);

        deltax_ = CalcXStep(flow_parameters_, imax_);
        deltay_ = CalcYStep(flow_parameters_, jmax_);
    }


    SupersonicCone::~SupersonicCone()
    {
    }

    double SupersonicCone::CalcXStep(const FlowParameters& params, int size) {
        return params.plate_length / size;
    }

    double SupersonicCone::CalcYStep(const FlowParameters& params, int size) {
        double rho = params.P_inf / (params.T_inf * params.R);
        double u = params.M_inf * params.a_inf;
        double Re = rho * u * params.plate_length / params.mu;
        double delta = 5 * params.plate_length / sqrt(Re);

        double lvert = 5 * delta;
        return lvert / size;
    }

    void SupersonicCone::Run() {

        // Sets initial conditions.
        InitializeFlowFliedVariables();

        int mod = 1;

        for (int it = 0; it < maxit; it++) {

            if (it % mod == 0)
                std::cout << "Iteration " << it << " ";

            double delta_t = CalcTStep(imax_, jmax_, deltax_, deltay_, flow_parameters_, u_, v_, rho_, P_, T_, 0.6);

            rho_old_ = rho_;
            u_old_ = u_;
            v_old_ = v_;
            P_old_ = P_;

            maccormack_solver_.UpdatePredictor(delta_t, deltax_, deltay_, imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

            BoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

            maccormack_solver_.UpdateCorrector(delta_t, deltax_, deltay_, imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

            BoundaryConditions(imax_, jmax_, flow_parameters_, u_, v_, rho_, P_, T_, e_);

            double diff = 0.0;
            if (CheckConvergence(diff)) {
                break;
            }

            if (it % mod == 0)
                std::cout << diff << std::endl;
        }

        WriteInFile(P_, deltax_, deltay_, "Pressure");
        WriteInFile(rho_, deltax_, deltay_, "Density");
        WriteInFile(T_, deltax_, deltay_, "Temperature");
        WriteInFile(u_, deltax_, deltay_, "VelocityX");
        WriteInFile(v_, deltax_, deltay_, "VelocityY");
    }

    void SupersonicCone::InitializeFlowFliedVariables() {

        for (int i = 0; i < imax_; i++) {
            for (int j = 0; j < jmax_; j++) {
                T_.Get(i, j) = flow_parameters_.T_inf;
                P_.Get(i, j) = flow_parameters_.P_inf;
                rho_.Get(i, j) = flow_parameters_.P_inf / flow_parameters_.T_inf / flow_parameters_.R;
                u_.Get(i, j) = flow_parameters_.M_inf * flow_parameters_.a_inf;
                v_.Get(i, j) = 0;
                M_.Get(i, j) = flow_parameters_.M_inf;
                e_.Get(i, j) = flow_parameters_.cv * flow_parameters_.T_inf;

                outside_.Set(i, j, !(i > margin && j < height));
            }
        }

        // Wall conditions

        /*
        for (int i = 20; i < imax_; i++) {
            for (int j = 0; j < jmax_; j++) {
                if (j == 0 ||atan2(j*deltay_, (i-20)*deltax_) < angle / 180 * pi) {
                    T_.Get(i, j) = flow_parameters_.T_wall;
                    u_.Get(i, j) = 0;
                    v_.Get(i, j) = 0;
                    M_.Get(i, j) = 0;
                    rho_.Get(i, j) = flow_parameters_.P_inf / T_.Get(i, 0) / flow_parameters_.R;
                    e_.Get(i, j) = flow_parameters_.cv * flow_parameters_.T_wall;
                }
            }
        }
        */

        for (int i = margin; i < imax_; i++) {
            T_.Get(i, height) = flow_parameters_.T_wall;
            u_.Get(i, height) = 0;
            v_.Get(i, height) = 0;
            M_.Get(i, height) = 0;
            rho_.Get(i, height) = flow_parameters_.P_inf / T_.Get(i, 0) / flow_parameters_.R;
            e_.Get(i, height) = flow_parameters_.cv * flow_parameters_.T_wall;
        }

        for (int j = 0; j <= height; j++) {
            T_.Get(margin, j) = flow_parameters_.T_wall;
            u_.Get(margin, j) = 0;
            v_.Get(margin, j) = 0;
            M_.Get(margin, j) = 0;
            rho_.Get(margin, j) = flow_parameters_.P_inf / flow_parameters_.T_wall / flow_parameters_.R;
            e_.Get(margin, j) = flow_parameters_.cv * flow_parameters_.T_wall;
        }

    }

    void SupersonicCone::BoundaryConditions(int imax, int jmax, const FlowParameters& params,
        Array2D<double>& u, Array2D<double>& v, Array2D<double>& rho, Array2D<double>& P,
        Array2D<double>& T, Array2D<double>& e) {

        for (int i = margin + 1; i < imax_; i++) {
            for (int j = 0; j < height; j++) {
                T_.Get(i, j) = flow_parameters_.T_inf;
                P_.Get(i, j) = flow_parameters_.P_inf;
                rho_.Get(i, j) = flow_parameters_.P_inf / flow_parameters_.T_inf / flow_parameters_.R;
                u_.Get(i, j) = flow_parameters_.M_inf * flow_parameters_.a_inf;
                v_.Get(i, j) = 0;
                M_.Get(i, j) = flow_parameters_.M_inf;
                e_.Get(i, j) = flow_parameters_.cv * flow_parameters_.T_inf;

                outside_.Set(i, j, !(i > margin && j < height));
            }
        }

        for (int j = 0; j < jmax; j++) {
            u.Set(0, j, params.a_inf*params.M_inf);
            v.Set(0, j, 0);
            P.Set(0, j, params.P_inf);
            T.Set(0, j, params.T_inf);
            rho.Set(0, j, params.P_inf / params.R / params.T_inf);
            e.Set(0, j, params.T_inf*params.cv);
        }

        for (int i = 0; i < imax; i++) {
            u.Set(i, jmax - 1, params.a_inf*params.M_inf);
            v.Set(i, jmax - 1, 0);
            P.Set(i, jmax - 1, params.P_inf);
            T.Set(i, jmax - 1, params.T_inf);
            rho.Set(i, jmax - 1, params.P_inf / params.R / params.T_inf);
            e.Set(i, jmax - 1, params.T_inf*params.cv);
        }

        double rho_ij;
        double u_ij;
        double p_ij;
        
        for (int i = 0; i < margin; i++) {
            rho_ij = 2 * rho.Get(i, 1) - rho.Get(i, 2);
            u_ij = 2 * u.Get(i, 1) - u.Get(i, 2);
            p_ij = 2 * P.Get(i, 1) - P.Get(i, 2);

            u.Set(i, 0, u_ij);
            v.Set(i, 0, 0);
            rho.Set(i, 0, rho_ij);
            P.Set(i, 0, p_ij);
            T.Set(i, 0, p_ij / params.R / rho_ij);
            e.Set(i, 0, T.Get(i, 0)*params.cv);
        }
        

        for (int j = 0; j <= height; j++) {
            p_ij = 2 * P.Get(margin + 1, j) - P.Get(margin + 1, j);
            T_.Get(margin, j) = flow_parameters_.T_wall;
            u_.Get(margin, j) = 0;
            v_.Get(margin, j) = 0;
            P.Set(margin, j, p_ij);
            rho_.Get(margin, j) = flow_parameters_.P_inf / flow_parameters_.T_wall / flow_parameters_.R;
            e_.Get(margin, j) = flow_parameters_.cv * flow_parameters_.T_wall;
        }
        /*
        for (int i = margin + 1; i < imax_; i++) {
            bool vis = false;
            for (int j = jmax_ - 1; j >= 0; j--) {
                if (j == 0 || atan2(j*deltay_, (i - 20)*deltax_) < 15 / 180 * pi) {
                    p_ij = 2 * P.Get(i, j + 1) - P.Get(i, j + 2);

                    u.Set(i, j, 0);
                    v.Set(i, j, 0);
                    P.Set(i, j, p_ij);
                    T.Set(i, j, params.T_wall);
                    rho.Set(i, j, p_ij / params.R / params.T_wall);
                    e.Set(i, j, params.T_wall*params.cv); 
                }
            }
        }
        */

        for (int j = 0; j < jmax; j++) {
            u.Set(imax - 1, j, 2 * u.Get(imax - 2, j) - u.Get(imax - 3, j));
            v.Set(imax - 1, j, 2 * v.Get(imax - 2, j) - v.Get(imax - 3, j));
            P.Set(imax - 1, j, 2 * P.Get(imax - 2, j) - P.Get(imax - 3, j));
            T.Set(imax - 1, j, 2 * T.Get(imax - 2, j) - T.Get(imax - 3, j));
            rho.Set(imax - 1, j, 2 * rho.Get(imax - 2, j) - rho.Get(imax - 3, j));
            e.Set(imax - 1, j, 2 * e.Get(imax - 2, j) - e.Get(imax - 3, j));;
        }







        for (int i = margin; i < imax_; i++) {
            p_ij = 2 * P.Get(i, height + 1) - P.Get(i, height + 2);
            T_.Get(i, height) = flow_parameters_.T_wall;
            u_.Get(i, height) = 0;
            v_.Get(i, height) = 0;
            P.Set(i, height, p_ij);
            rho_.Get(i, height) = flow_parameters_.P_inf / flow_parameters_.T_wall / flow_parameters_.R;
            e_.Get(i, height) = flow_parameters_.cv * flow_parameters_.T_wall;
        }

    }

    bool SupersonicCone::CheckConvergence(double& diff) {

        for (int i = 0; i < imax_; i++) {
            for (int j = 0; j < jmax_; j++) {
                double tmp = abs(rho_.Get(i, j) - rho_old_.Get(i, j));
                tmp = std::max(abs(u_.Get(i, j) - u_old_.Get(i, j)), tmp);
                tmp = std::max(abs(v_.Get(i, j) - v_old_.Get(i, j)), tmp);
                tmp = std::max(abs(P_.Get(i, j) - P_old_.Get(i, j)), tmp);
                diff = std::max(diff, tmp);
            }
        }

        return diff < 1.0e-6;
    }
}