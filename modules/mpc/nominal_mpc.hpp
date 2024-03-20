#ifndef _NOMINAL_MPC_HPP
#define _NOMINAL_MPC_HPP

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>
#include <Eigen/Dense>
#include <chrono>

#include "quadrotor_dynamics/nominal_quad_dynamic.hpp"

using namespace std;
using namespace Eigen;

int cnt = 0;

template <int _dt_num, int _dt_den, int _n_step>
class NominalMpc_ : protected NominalQuadDynamic {
public:
    static constexpr int x_dim_ = 10;
    static constexpr int u_dim_ = 4;

private:
    Matrix<double, x_dim_, 1> state_[_n_step];
    Matrix<double, x_dim_, _n_step * u_dim_> state_g_[_n_step];
    Matrix<double, _n_step, x_dim_> ref_state_;
    Matrix<double, _n_step, u_dim_> ref_u_;
    Matrix<double, x_dim_, 1> init_state_;
    Matrix<double, x_dim_ + u_dim_, 1> cost_w_;
    const double dt_;
    nlopt::opt opt;

public:
    NominalMpc_(double hover_ratio, string algorithm = "LD_LBFGS", int maxeval = 50) : 
        NominalQuadDynamic(hover_ratio),
        dt_((double)_dt_num / _dt_den),
        opt(algorithm.c_str(), u_dim_ * _n_step) {
        for (int k = 0; k < _n_step; k++) {
            state_[k].setZero();
            state_g_[k].setZero();
        }
        vector<double> lb(u_dim_ * _n_step), ub(u_dim_ * _n_step);
        //boundaries on optimization variables
        for (int k = 0; k < _n_step; k++) {
            lb[k * u_dim_ + 0] = -M_PI;
            lb[k * u_dim_ + 1] = -M_PI;
            lb[k * u_dim_ + 2] = -M_PI;
            lb[k * u_dim_ + 3] = 0;
            ub[k * u_dim_ + 0] = M_PI;
            ub[k * u_dim_ + 1] = M_PI;
            ub[k * u_dim_ + 2] = M_PI;
            ub[k * u_dim_ + 3] = 1;
        }
        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);
        opt.set_min_objective((nlopt::vfunc)&NominalMpc_::cost_func, this);
        opt.set_xtol_rel(1e-24);
        opt.set_ftol_rel(1e-12);
        opt.set_maxeval(maxeval);
        // opt.set_vector_storage(1);
    }

    static double cost_func(const vector<double> &u, vector<double> &grad, NominalMpc_ *instance)
    {
        Matrix<double, x_dim_, 1> *state = instance->state_;
        Matrix<double, x_dim_, _n_step * u_dim_> *state_g = instance->state_g_;
        const double &dt = instance->dt_;
        VectorXd g = VectorXd::Zero(u.size());
        VectorXd uvec(u_dim_);

        //update preditive model
        uvec << u[0], u[1], u[2], u[3];
        Matrix<double, x_dim_, x_dim_> x1dotx0;
        Matrix<double, x_dim_, u_dim_> x1dotu;
        x1dotx0.setZero();
        x1dotu.setZero();
        instance->rk4_func(instance->init_state_, uvec, dt, state[0], x1dotx0, x1dotu);
        state_g[0].block(0, 0, x_dim_, u_dim_) = x1dotu;
        for (int k = 1; k < _n_step; k++) {
            uvec << u[0 + k * u_dim_], u[1 + k * u_dim_], u[2 + k * u_dim_], u[3 + k * u_dim_];
            x1dotx0.setZero();
            x1dotu.setZero();
            instance->rk4_func(state[k - 1], uvec, dt, state[k], x1dotx0, x1dotu);
            state_g[k].block(0, k * u_dim_, x_dim_, u_dim_) = x1dotu;
            state_g[k].block(0, 0, x_dim_, u_dim_ * k) = x1dotx0 * state_g[k - 1].block(0, 0, x_dim_, u_dim_ * k);
        }

        //calculate cost
        Matrix<double, _n_step, x_dim_> &ref_state = instance->ref_state_;
        Matrix<double, _n_step, u_dim_> &ref_u = instance->ref_u_;
        Matrix<double, x_dim_ + u_dim_, 1> &cost_w = instance->cost_w_;
        double pvcost = 0.0;
        for (int k = 0; k < _n_step; k++) {
            double scaling = 1.0;
            if (k < _n_step - 1) {
                scaling = 0.1;
            } else {
                scaling = 1.0;
            }
            pvcost += scaling * (cost_w[0] * pow(state[k](0) - ref_state(k, 0), 2)
                + cost_w[1] * pow(state[k](1) - ref_state(k, 1), 2)
                + cost_w[2] * pow(state[k](2) - ref_state(k, 2), 2)
                + cost_w[3] * pow(state[k](3) - ref_state(k, 3), 2)
                + cost_w[4] * pow(state[k](4) - ref_state(k, 4), 2)
                + cost_w[5] * pow(state[k](5) - ref_state(k, 5), 2));
            g += scaling * 2 * (
                  cost_w[0] * (state[k](0) - ref_state(k, 0)) * state_g[k].row(0).transpose()
                + cost_w[1] * (state[k](1) - ref_state(k, 1)) * state_g[k].row(1).transpose()
                + cost_w[2] * (state[k](2) - ref_state(k, 2)) * state_g[k].row(2).transpose()
                + cost_w[3] * (state[k](3) - ref_state(k, 3)) * state_g[k].row(3).transpose()
                + cost_w[4] * (state[k](4) - ref_state(k, 4)) * state_g[k].row(4).transpose()
                + cost_w[5] * (state[k](5) - ref_state(k, 5)) * state_g[k].row(5).transpose());
        }
        double ucost = 0.0;
        for (int k = 0; k < _n_step; k++) {
            double scaling = 1.0;
            if (k < _n_step - 1) {
                scaling = 0.1;
            } else {
                scaling = 0.1;
            }
            int of = k * u_dim_;
            ucost += scaling * (cost_w[10] * pow(u[of + 0] - ref_u(k, 0), 2)
                + cost_w[11] * pow(u[of + 1] - ref_u(k, 1), 2)
                + cost_w[12] * pow(u[of + 2] - ref_u(k, 2), 2)
                + cost_w[13] * pow(u[of + 3] - ref_u(k, 3), 2));
            g[of + 0] += scaling * 2 * cost_w[10] * (u[of + 0] - ref_u(k, 0));
            g[of + 1] += scaling * 2 * cost_w[11] * (u[of + 1] - ref_u(k, 1));
            g[of + 2] += scaling * 2 * cost_w[12] * (u[of + 2] - ref_u(k, 2));
            g[of + 3] += scaling * 2 * cost_w[13] * (u[of + 3] - ref_u(k, 3));
        }

        for (int i = 0; i < grad.size(); i++) {
            grad[i] = g[i];
        }

        return pvcost + ucost;
    }

    void set_w(Matrix<double, x_dim_ + u_dim_, 1> &w) {
        cost_w_ = w;
    }

    int solve(const Matrix<double, x_dim_, 1> &state,
        const Matrix<double, _n_step, x_dim_> &ref_state,
        Matrix<double, _n_step, u_dim_> &u, 
        Matrix<double, _n_step, x_dim_> &x_predict,
        double &solve_time) {
        vector<double> uv(u_dim_ * _n_step);
        //initial solution
        for (int k = 0; k < _n_step; k++) {
            uv[k * u_dim_ + 0] = u(k, 0);
            uv[k * u_dim_ + 1] = u(k, 1);
            uv[k * u_dim_ + 2] = u(k, 2);
            uv[k * u_dim_ + 3] = u(k, 3);
        }
        init_state_ = state;
        ref_state_ = ref_state;
        for (int k = 0; k < _n_step; k++) {
            ref_u_.row(k) << 0, 0, 0, this->hover_ratio_;
        }

        double minf;
        try{
            auto beforeTime = chrono::steady_clock::now();
            opt.optimize(uv, minf);
            auto afterTime = chrono::steady_clock::now();
            solve_time = chrono::duration<double>(afterTime - beforeTime).count();

            for (int k = 0; k < _n_step; k++) {
                u(k, 0) = uv[k * u_dim_ + 0];
                u(k, 1) = uv[k * u_dim_ + 1];
                u(k, 2) = uv[k * u_dim_ + 2];
                u(k, 3) = uv[k * u_dim_ + 3];
            }
            for (int k = 0; k < _n_step; k++) {
                x_predict.row(k) = state_[k];
            }
            // cout << "Solve time: " << solve_time * 1e6 << " us" << endl;

            return EXIT_SUCCESS;
        } catch(exception &e) {
            cerr << "nlopt failed: " << e.what() << endl;
            return EXIT_FAILURE;
        }
    }
};

#endif