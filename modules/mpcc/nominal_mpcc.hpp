#ifndef _NOMINAL_MPCC_HPP
#define _NOMINAL_MPCC_HPP

#define USE_EXTENDED_DYNAMICS 0

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>
#include <Eigen/Dense>
#include <chrono>
#include <list>

#if USE_EXTENDED_DYNAMICS
#include "quadrotor_dynamics/extended_quad_dynamic.hpp"
#else
#include "quadrotor_dynamics/nominal_quad_dynamic.hpp"
#endif
#include "bspline/uniform_bspline.hpp"
#include "map/sdf.hpp"

using namespace std;
using namespace Eigen;

#if USE_EXTENDED_DYNAMICS
class NominalMpcc : protected ExtendedQuadDynamic {
#else
class NominalMpcc : protected NominalQuadDynamic {
#endif
public:
#if USE_EXTENDED_DYNAMICS
    static constexpr int x_dim_ = ExtendedQuadDynamic::x_dim_;
#else
    static constexpr int x_dim_ = NominalQuadDynamic::x_dim_;
#endif
    static constexpr int u_dim_ = 4;
    static constexpr int _dt_num = 1;
    static constexpr int _dt_den = 10;
    static constexpr int _n_step = 10;

private:
    Matrix<double, x_dim_, 1> state_[_n_step];
    Matrix<double, x_dim_, _n_step * u_dim_> state_g_[_n_step];
    Matrix<double, 3, 1> acc_[_n_step];
    Matrix<double, 3, _n_step * u_dim_> acc_g_[_n_step];
    // Matrix<double, _n_step, x_dim_> ref_state_;
    const SdfMap *sdf_;
    MatrixXd ctrl_pts_;
    MatrixXd v_ctrl_pts_;
    MatrixXd a_ctrl_pts_;
    double ts_;
    double t_max_;
    Matrix<double, u_dim_, 1> last_u_;
    Matrix<double, _n_step, u_dim_> ref_u_;
    Matrix<double, x_dim_, 1> init_state_;
    Matrix<double, 3 + u_dim_ + 1 + u_dim_ + 1, 1> cost_w; //contour_err_w, lag_err_w, time_w, u_w
    vector<double> tmp_u_;
    bool flag;
    const double dt_;
    Vector3d disturbance_acc_;
    const vector<DynObs> *dynobs_;

    list<Vector3d> past_disturbances_;

    double cbf_cost_;
    double t_w_ratio_;

    nlopt::opt opt_;

public:
    NominalMpcc(double hover_ratio, string algorithm = "LD_LBFGS", int maxeval = 50);

    static void v_b_constraint(unsigned m, double *result, unsigned n, const double *u,
                             double *gradient, /* NULL if not needed */
                             NominalMpcc *instance);

    static void collision_constraint(unsigned m, double *result, unsigned n, const double *u,
                             double *gradient, /* NULL if not needed */
                             NominalMpcc *instance);

    static void cbf_constraint(unsigned m, double *result, unsigned n, const double *u,
                             double *gradient, /* NULL if not needed */
                             NominalMpcc *instance);

    static double cost_func(const vector<double> &u, vector<double> &grad, NominalMpcc *instance);

    void set_w(Matrix<double, 3 + u_dim_ + 1 + u_dim_ + 1, 1> &w);

    int solve(const Matrix<double, x_dim_, 1> &state,
        const SdfMap &sdf,
        const MatrixXd &ctrl_pts,
        const MatrixXd &v_ctrl_pts,
        const MatrixXd &a_ctrl_pts,
        const double ts,
        const double len,
        const Matrix<double, u_dim_, 1> last_u,
        const Vector3d disturbance_acc,
        const vector<DynObs> &dynobs,
        Matrix<double, _n_step, u_dim_> &u, 
        Matrix<double, _n_step, x_dim_> &x_predict,
        Matrix<double, _n_step + 1, 1> &t_index,
        double &solve_time);
};

#endif
