#include "nominal_mpcc.hpp"

int cnt = 0;

VectorXd vector2Vector(const vector<double> &v) {
    VectorXd vv(v.size());
    for(int i = 0; i < v.size(); i++) {
        vv[i] = v[i];
    }
    return vv;
}

NominalMpcc::NominalMpcc(double hover_ratio, string algorithm, int maxeval) : 
#if USE_EXTENDED_DYNAMICS
    ExtendedQuadDynamic(hover_ratio),
#else
    NominalQuadDynamic(hover_ratio),
#endif
    dt_((double)_dt_num / _dt_den),
    opt_(algorithm.c_str(), u_dim_ * _n_step + _n_step + 1) {
    for (int k = 0; k < _n_step; k++) {
        state_[k].setZero();
        state_g_[k].setZero();
        acc_[k].setZero();
        acc_g_[k].setZero();
    }
    flag = false;
    opt_.set_min_objective((nlopt::vfunc)&NominalMpcc::cost_func, this);
    vector<double> tols(_n_step);
    for (int i = 0; i < tols.size(); i++) {
        tols[i] = 1e-8;
    }
    opt_.add_inequality_mconstraint((nlopt::mfunc)&NominalMpcc::v_b_constraint, this, tols);
    opt_.set_xtol_rel(1e-16);
    opt_.set_ftol_rel(1e-10);
    opt_.set_maxeval(maxeval);
    nlopt::opt local_opt("LD_LBFGS", u_dim_ * _n_step + _n_step + 1);
    local_opt.set_xtol_rel(1e-16);
    local_opt.set_ftol_rel(1e-10);
    local_opt.set_maxeval(maxeval);
    // local_opt.set_vector_storage(16);
    opt_.set_local_optimizer(local_opt);
}

//velocity constraints
void NominalMpcc::v_b_constraint(unsigned m, double *result, unsigned n, const double *u,
                            double *gradient, /* NULL if not needed */
                            NominalMpcc *instance) {
    Matrix<double, x_dim_, 1> *state = instance->state_;
    Matrix<double, x_dim_, _n_step * u_dim_> *state_g = instance->state_g_;
    const double &dt = instance->dt_;

    for (int k = 0; k < m; k++) {
        result[k] = pow(state[k][3], 2) + pow(state[k][4], 2) + pow(state[k][5], 2) - pow(20.0, 2);
        if (gradient) {
            for (int i = 0; i < _n_step * u_dim_; i++) {
                gradient[k * n + i] = 2 * state[k][3] * state_g[k](3, i) + 
                                        2 * state[k][4] * state_g[k](4, i) +
                                        2 * state[k][5] * state_g[k](5, i);
            }
            for (int i = _n_step * u_dim_; i < n; i++) {
                gradient[k * n + i] = 0;
            }
        }
    }
}

double NominalMpcc::cost_func(const vector<double> &u, vector<double> &grad, NominalMpcc *instance)
{
    Matrix<double, x_dim_, 1> *state = instance->state_;
    Matrix<double, x_dim_, _n_step * u_dim_> *state_g = instance->state_g_;
    Matrix<double, 3, 1> *acc = instance->acc_;
    Matrix<double, 3, _n_step * u_dim_> *acc_g = instance->acc_g_;
    const double &dt = instance->dt_;
    Matrix<double, _n_step, u_dim_> &ref_u = instance->ref_u_;
    auto &cost_w = instance->cost_w;

    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> g;
    g.setZero();
    Matrix<double, u_dim_, 1> uvec;

    //update preditive model using RK4 
    //Compared with Euler method, RK4 has better accuracy in solving the state equation
    uvec << u[0], u[1], u[2], u[3];
    Matrix<double, x_dim_, x_dim_> x1dotx0;
    Matrix<double, x_dim_, u_dim_> x1dotu;
    Matrix<double, 3, x_dim_> accdotx0;
    Matrix<double, 3, u_dim_> accdotu;
    instance->rk4_func(instance->init_state_, uvec, instance->disturbance_acc_, dt, state[0], x1dotx0, x1dotu, acc[0], accdotx0, accdotu);
    state_g[0].block(0, 0, x_dim_, u_dim_) = x1dotu;
    acc_g[0].block(0, 0, 3, u_dim_) = accdotu;
    for (int k = 1; k < _n_step; k++) {
        uvec << u[0 + k * u_dim_], u[1 + k * u_dim_], u[2 + k * u_dim_], u[3 + k * u_dim_];
        instance->rk4_func(state[k - 1], uvec, k < 2 ? instance->disturbance_acc_ : Vector3d(0, 0, 0), dt, state[k], x1dotx0, x1dotu, acc[k], accdotx0, accdotu);
        state_g[k].block(0, k * u_dim_, x_dim_, u_dim_) = x1dotu;
        state_g[k].block(0, 0, x_dim_, u_dim_ * k) = x1dotx0 * state_g[k - 1].block(0, 0, x_dim_, u_dim_ * k);
        acc_g[k].block(0, k * u_dim_, 3, u_dim_) = accdotu;
        acc_g[k].block(0, 0, 3, u_dim_ * k) = accdotx0 * state_g[k - 1].block(0, 0, x_dim_, u_dim_ * k);
    }

    //calculate reference state and tangent vector
    Matrix<double, 3, 1> ref_pos[_n_step + 1];
    Matrix<double, 3, _n_step + 1> ref_pos_g[_n_step + 1];
    Matrix<double, 3, 1> tangent[_n_step + 1];
    Matrix<double, 3, _n_step + 1> tangent_g[_n_step + 1];
    double t_sum = u[u_dim_ * _n_step + _n_step];
    for (int k = 0; k < _n_step + 1; k++) {
        bool over = false;
        if (t_sum >= instance->t_max_) {
            t_sum = instance->t_max_;
            over = true;
        }
        Matrix<double, 3, 1> ref_pos_k_g;
        ref_pos[k] = UniformBspline::getBsplineValueFast(instance->ts_, instance->ctrl_pts_, t_sum, 3, &ref_pos_k_g);
        Matrix<double, 3, 1> tmp2;
        Matrix<double, 3, 1> tmp = UniformBspline::getBsplineValueFast(instance->ts_, instance->v_ctrl_pts_, t_sum - instance->ts_, 2, &tmp2);
        double snorm = tmp.squaredNorm();
        double norm = tmp.norm();
        tangent[k] = tmp / norm;
        ref_pos_g[k].setZero();
        tangent_g[k].setZero();
        if (!over) {
            ref_pos_g[k].block(0, _n_step, 3, 1) = ref_pos_k_g;
            tangent_g[k].block(0, _n_step, 3, 1) = (tmp2 * norm - tmp * 0.5 * (1 / norm) * 2 * tmp.dot(tmp2)) / snorm;
            for (int i = 0; i < k; i++) {
                ref_pos_g[k].block(0, i, 3, 1) = ref_pos_k_g;
                tangent_g[k].block(0, i, 3, 1) = (tmp2 * norm - tmp * 0.5 * (1 / norm) * 2 * tmp.dot(tmp2)) / snorm;
            }
        }
        if (k < _n_step) {
            t_sum += u[u_dim_ * _n_step + k];
        }
    }

    //calculate error vector between the predicted state and the reference state
    Matrix<double, 3, 1> err[_n_step + 1];
    Matrix<double, 3, u_dim_ * _n_step + _n_step + 1> err_g[_n_step + 1];
    for (int k = 0; k < _n_step + 1; k++) {
        if (k > 0) {
            err[k] = state[k - 1].block(0, 0, 3, 1) - ref_pos[k];
            err_g[k].setZero();
            err_g[k].block(0, u_dim_ * _n_step, 3, _n_step + 1) = -ref_pos_g[k];
            err_g[k].block(0, 0, 3, u_dim_ * _n_step) = state_g[k - 1].block(0, 0, 3, u_dim_ * _n_step);
        } else {
            err[k] = instance->init_state_.block(0, 0, 3, 1) - ref_pos[k];
            err_g[k].setZero();
            err_g[k].block(0, u_dim_ * _n_step, 3, _n_step + 1) = -ref_pos_g[k];
        }
    }

    //calculate lag error
    Matrix<double, 3, 1> lag_err[_n_step + 1];
    Matrix<double, 3, u_dim_ * _n_step + _n_step + 1> lag_err_g[_n_step + 1];
    double lag_e[_n_step + 1];
    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> lag_e_g[_n_step + 1];
    for (int k = 0; k < _n_step + 1; k++) {
        lag_e[k] = err[k].dot(tangent[k]);
        lag_e_g[k] = (tangent[k].transpose() * err_g[k]).transpose();
        lag_e_g[k].block(u_dim_ * _n_step, 0, _n_step + 1, 1) += (err[k].transpose() * tangent_g[k]).transpose();
        lag_err[k] = lag_e[k] * tangent[k];
        lag_err_g[k] = tangent[k] * lag_e_g[k].transpose();
        lag_err_g[k].block(0, u_dim_ * _n_step, 3, _n_step + 1) += lag_e[k] * tangent_g[k];
        lag_e_g[k] = 2 * lag_e[k] * lag_e_g[k];
        lag_e[k] = pow(lag_e[k], 2);
    }
    
    //calculate contouring error
    Matrix<double, 3, 1> contour_err[_n_step + 1];
    Matrix<double, 3, u_dim_ * _n_step + _n_step + 1> contour_err_g[_n_step + 1];
    double contour_e[_n_step + 1];
    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> contour_e_g[_n_step + 1];
    for (int k = 0; k < _n_step + 1; k++) {
        contour_err[k] = err[k] - lag_err[k];
        contour_err_g[k] = err_g[k] - lag_err_g[k];
        contour_e[k] = contour_err[k].squaredNorm();
        contour_e_g[k] = (2 * contour_err[k].transpose() * contour_err_g[k]).transpose();
    }

    //calculate the cost of tracking progress
    double t_cost = -u[u_dim_ * _n_step + _n_step];
    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> t_cost_g;
    t_cost_g.setZero();
    t_cost_g[u_dim_ * _n_step + _n_step] = -1;
    for (int k = 0; k < _n_step; k++) {
        t_cost -= u[u_dim_ * _n_step + k];
        t_cost_g[u_dim_ * _n_step + k] = -1;
    }
    // t_cost -= u[u_dim_ * _n_step] * 0.5;
    // t_cost_g[u_dim_ * _n_step] = -1 * 0.5;

    //calculate the cost of violating distance constraints (set its weight to 0 when using CBF constraints)
    double c_cost = 0.0;
    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> c_cost_g;
    c_cost_g.setZero();
    for (int k = 0; k < _n_step; k++) {
        auto dis = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k].block(0, 0, 3, 1));
        if (dis.first < 0.25) {
            double tmp = dis.first - 0.25;
            double tmp2 = pow(tmp, 2);
            c_cost += tmp2;
            c_cost_g.block(0, 0, u_dim_ * _n_step, 1) += (2 * tmp * dis.second.transpose() * state_g[k].block(0, 0, 3, u_dim_ * _n_step)).transpose();
        }
    }
    // cout << endl;

    //calculate the cost of violating CBF constraints
    double cbf_cost = 0.;
    // double gamma = 0.7;
    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> cbf_cost_g;
    cbf_cost_g.setZero();
    bool cbf_en = true;
#if 1
    for (int k = 0; k < _n_step; k++) {
        std::pair<double, Eigen::Vector3d> dis1;
        if (k == 0) {
            dis1 = instance->sdf_->get_dist_with_grad_trilinear<double>(instance->init_state_.block(0, 0, 3, 1));
        } else {
            dis1 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k - 1].block(0, 0, 3, 1));
        }
        auto dis2 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k].block(0, 0, 3, 1));
        double gamma = cbf_en ? 0.65 : 1.0;//0.2 * k / (_n_step - 1) + 0.8;
        double cbf = (1 - gamma) * (dis1.first - 0.25) - (dis2.first - 0.25);
        if (cbf > 0.0) {
            double tmp = cbf;//pow(cbf, 2);
            double scale = cbf_en ? pow(0.9, k) : 1.0;
            cbf_cost += scale * tmp;
            if (k == 0) {
                cbf_cost_g.block(0, 0, u_dim_ * _n_step, 1) -= scale * (dis2.second.transpose() * state_g[k].block(0, 0, 3, u_dim_ * _n_step)).transpose();
            } else {
                cbf_cost_g.block(0, 0, u_dim_ * _n_step, 1) += scale * (((1 - gamma) * dis1.second.transpose() * state_g[k - 1].block(0, 0, 3, u_dim_ * _n_step)
                    - dis2.second.transpose() * state_g[k].block(0, 0, 3, u_dim_ * _n_step))).transpose();
            }
        }
    }
#else
    for (int k = 0; k < _n_step; k++) {
        double gamma1 = cbf_en ? 0.65 : 1.0;
        double gamma2 = cbf_en ? 0.9 : 1.0;
        double gamma3 = cbf_en ? 0.95 : 1.0;
        double scale = cbf_en ? pow(0.9, k) : 1.0;
        if (k == 0) {
            std::pair<double, Eigen::Vector3d> dis1 = instance->sdf_->get_dist_with_grad_trilinear<double>(instance->init_state_.block(0, 0, 3, 1));
            auto dis2 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k].block(0, 0, 3, 1));
            double cbf = (dis2.first - 0.25) + (gamma1 - 1) * (dis1.first - 0.25);
            if (cbf < 0.0) {
                double tmp = -cbf;
                cbf_cost += scale * tmp;
                cbf_cost_g.block(0, 0, u_dim_ * _n_step, 1) -= scale * (dis2.second.transpose() * state_g[k].block(0, 0, 3, u_dim_ * _n_step)).transpose();
            }
        } else if (k == 1) {
            std::pair<double, Eigen::Vector3d> dis1 = instance->sdf_->get_dist_with_grad_trilinear<double>(instance->init_state_.block(0, 0, 3, 1));
            std::pair<double, Eigen::Vector3d> dis2 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k - 1].block(0, 0, 3, 1));
            auto dis3 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k].block(0, 0, 3, 1));
            double cbf = (dis3.first - 0.25) 
                        + (gamma1 + gamma2 - 2) * (dis2.first - 0.25) 
                        + (gamma1 - 1) * (gamma2 - 1) * (dis1.first - 0.25);
            if (cbf < 0.0) {
                double tmp = -cbf;
                cbf_cost += scale * tmp;
                cbf_cost_g.block(0, 0, u_dim_ * _n_step, 1) -= scale * (dis3.second.transpose() * state_g[k].block(0, 0, 3, u_dim_ * _n_step)
                                                            + (gamma1 + gamma2 - 2) * dis2.second.transpose() * state_g[k - 1].block(0, 0, 3, u_dim_ * _n_step)).transpose();
            }
        } else if (k == 2) {
            std::pair<double, Eigen::Vector3d> dis1 = instance->sdf_->get_dist_with_grad_trilinear<double>(instance->init_state_.block(0, 0, 3, 1));
            std::pair<double, Eigen::Vector3d> dis2 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k - 2].block(0, 0, 3, 1));
            std::pair<double, Eigen::Vector3d> dis3 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k - 1].block(0, 0, 3, 1));
            auto dis4 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k].block(0, 0, 3, 1));
            double cbf = (dis4.first - 0.25) 
                        + (gamma1 + gamma2 + gamma3 - 3) * (dis3.first - 0.25) 
                        + ((gamma1 - 1) * (gamma2 + gamma3 - 2) + (gamma2 - 1) * (gamma3 - 1)) * (dis2.first - 0.25)
                        + (gamma1 - 1) * (gamma2 - 1) * (gamma3 - 1) * (dis1.first - 0.25);
            if (cbf < 0.0) {
                double tmp = -cbf;
                cbf_cost += scale * tmp;
                cbf_cost_g.block(0, 0, u_dim_ * _n_step, 1) -= scale * (
                    dis4.second.transpose() * state_g[k].block(0, 0, 3, u_dim_ * _n_step)
                    + (gamma1 + gamma2 + gamma3 - 3) * dis3.second.transpose() * state_g[k - 1].block(0, 0, 3, u_dim_ * _n_step)
                    + ((gamma1 - 1) * (gamma2 + gamma3 - 2) + (gamma2 - 1) * (gamma3 - 1)) * dis2.second.transpose() * state_g[k - 2].block(0, 0, 3, u_dim_ * _n_step)).transpose();
            }
        } else {
            std::pair<double, Eigen::Vector3d> dis1 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k - 3].block(0, 0, 3, 1));
            std::pair<double, Eigen::Vector3d> dis2 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k - 2].block(0, 0, 3, 1));
            std::pair<double, Eigen::Vector3d> dis3 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k - 1].block(0, 0, 3, 1));
            auto dis4 = instance->sdf_->get_dist_with_grad_trilinear<double>(state[k].block(0, 0, 3, 1));
            double cbf = (dis4.first - 0.25) 
                        + (gamma1 + gamma2 + gamma3 - 3) * (dis3.first - 0.25) 
                        + ((gamma1 - 1) * (gamma2 + gamma3 - 2) + (gamma2 - 1) * (gamma3 - 1)) * (dis2.first - 0.25)
                        + (gamma1 - 1) * (gamma2 - 1) * (gamma3 - 1) * (dis1.first - 0.25);
            if (cbf < 0.0) {
                double tmp = -cbf;
                cbf_cost += scale * tmp;
                cbf_cost_g.block(0, 0, u_dim_ * _n_step, 1) -= scale * (
                    dis4.second.transpose() * state_g[k].block(0, 0, 3, u_dim_ * _n_step)
                    + (gamma1 + gamma2 + gamma3 - 3) * dis3.second.transpose() * state_g[k - 1].block(0, 0, 3, u_dim_ * _n_step)
                    + ((gamma1 - 1) * (gamma2 + gamma3 - 2) + (gamma2 - 1) * (gamma3 - 1)) * dis2.second.transpose() * state_g[k - 2].block(0, 0, 3, u_dim_ * _n_step)
                    + (gamma1 - 1) * (gamma2 - 1) * (gamma3 - 1) * dis1.second.transpose() * state_g[k - 3].block(0, 0, 3, u_dim_ * _n_step)).transpose();
            }
        }
    }
#endif
    instance->cbf_cost_ = cbf_cost;

    //calculate heading angle cost
    double yawcost = 0.0;
    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> yawcost_g;
    yawcost_g.setZero();
    for (int k = 0; k < _n_step; k++) {
        double &vx = state[k][3];
        double &vy = state[k][4];
        double &qw = state[k][6];
        double &qx = state[k][7];
        double &qy = state[k][8];
        double &qz = state[k][9];
        double tmpx = 1. - 2. * (qy * qy + qz * qz);
        double tmpy = 2. * (qw * qz + qx * qy);
        double tmpnorm = pow(tmpx, 2) + pow(tmpy, 2);
        double yaw = std::atan2(tmpy, tmpx);
        Vector2d yaw_g_tmp(-tmpy / tmpnorm, tmpx / tmpnorm);
        Vector4d yaw_g;
        yaw_g.setZero();
        yaw_g(0) = yaw_g_tmp.y() * 2 * qz;
        yaw_g(1) = yaw_g_tmp.y() * 2 * qy;
        yaw_g(2) = yaw_g_tmp.x() * (-4 * qy) + yaw_g_tmp.y() * 2 * qx;
        yaw_g(3) = yaw_g_tmp.x() * (-4 * qz) + yaw_g_tmp.y() * 2 * qw;
        if (pow(vx, 2) + pow(vy, 2) > 1e-6) {
            double vang = std::atan2(vy, vx);
            double vxynorm = pow(vx, 2) + pow(vy, 2);
            Vector2d vang_g(-vy / vxynorm, vx / vxynorm);
            if (yaw - vang > M_PI) {
                yawcost += 2 * M_PI - (yaw - vang);
                yawcost_g.block(0, 0, u_dim_ * _n_step, 1) += 
                    (-yaw_g.transpose() * state_g[k].block(6, 0, 4, u_dim_ * _n_step)).transpose()
                    + (vang_g.transpose() * state_g[k].block(3, 0, 2, u_dim_ * _n_step)).transpose();
            } else if (vang - yaw > M_PI) {
                yawcost += 2 * M_PI + (yaw - vang);
                yawcost_g.block(0, 0, u_dim_ * _n_step, 1) += 
                    (yaw_g.transpose() * state_g[k].block(6, 0, 4, u_dim_ * _n_step)).transpose()
                    + (-vang_g.transpose() * state_g[k].block(3, 0, 2, u_dim_ * _n_step)).transpose();
            } else if (yaw > vang) {
                yawcost += yaw - vang;
                yawcost_g.block(0, 0, u_dim_ * _n_step, 1) += 
                    (yaw_g.transpose() * state_g[k].block(6, 0, 4, u_dim_ * _n_step)).transpose()
                    + (-vang_g.transpose() * state_g[k].block(3, 0, 2, u_dim_ * _n_step)).transpose();
            } else {
                yawcost += -yaw + vang;
                yawcost_g.block(0, 0, u_dim_ * _n_step, 1) += 
                    (-yaw_g.transpose() * state_g[k].block(6, 0, 4, u_dim_ * _n_step)).transpose()
                    + (vang_g.transpose() * state_g[k].block(3, 0, 2, u_dim_ * _n_step)).transpose();
            }
        }
    }

    //calculate the cost of violating velocity constraints
    double vcost = 0.0;
    Matrix<double, u_dim_ * _n_step + _n_step + 1, 1> vcost_g;
    vcost_g.setZero();
    for (int k = 0; k < _n_step; k++) {
        double tmp = pow(state[k][3], 2) + pow(state[k][4], 2) + pow(state[k][5], 2) - pow(15, 2);
        if (tmp > 0) {
            vcost += tmp;
            vcost_g.block(0, 0, u_dim_ * _n_step, 1) += 2 * state[k][3] * state_g[k].block(3, 0, 1, u_dim_ * _n_step).transpose()
                                                    + 2 * state[k][4] * state_g[k].block(4, 0, 1, u_dim_ * _n_step).transpose()
                                                    + 2 * state[k][5] * state_g[k].block(5, 0, 1, u_dim_ * _n_step).transpose();
        }
    }

    //calculate the cost of control input
    double ucost = 0.0;
    for (int k = 0; k < _n_step; k++) {
        int of = k * u_dim_;
        Matrix<double, u_dim_, 1> past_u;
        if (k == 0) {
            past_u = instance->last_u_;
        } else {
            past_u = Matrix<double, u_dim_, 1>(
                u[of - 4], u[of - 3], u[of - 2], u[of - 1]
            );
        }
        ucost += 
            // + cost_w[3] * pow(u[of + 0] - ref_u(k, 0), 2)
            // + cost_w[4] * pow(u[of + 1] - ref_u(k, 1), 2)
            + cost_w[6] * pow(u[of + 2] - ref_u(k, 2), 2)
            // + cost_w[6] * pow(u[of + 3] - ref_u(k, 3), 2)
            + cost_w[3] * pow(acc[k][0], 2)
            + cost_w[4] * pow(acc[k][1], 2)
            + cost_w[5] * pow(acc[k][2], 2)
            + cost_w[8] * pow(u[of + 0] - past_u(0), 2)
            + cost_w[9] * pow(u[of + 1] - past_u(1), 2)
            + cost_w[10] * pow(u[of + 2] - past_u(2), 2)
#if USE_EXTENDED_DYNAMICS
            + cost_w[11] * pow(u[of + 3], 2);
#else
            + cost_w[11] * pow(u[of + 3] - past_u(3), 2);
#endif
        // g[of + 0] += 2 * cost_w[3] * (u[of + 0] - ref_u(k, 0));
        // g[of + 1] += 2 * cost_w[4] * (u[of + 1] - ref_u(k, 1));
        g[of + 2] += 2 * cost_w[6] * (u[of + 2] - ref_u(k, 2));
        // g[of + 3] += 2 * cost_w[6] * (u[of + 3] - ref_u(k, 3));
        g.block(0, 0, u_dim_ * _n_step, 1) += 2 * cost_w[3] * acc[k][0] * acc_g[k].row(0).transpose();
        g.block(0, 0, u_dim_ * _n_step, 1) += 2 * cost_w[4] * acc[k][1] * acc_g[k].row(1).transpose();
        g.block(0, 0, u_dim_ * _n_step, 1) += 2 * cost_w[5] * acc[k][2] * acc_g[k].row(2).transpose();
        g[of + 0] += 2 * cost_w[8] * (u[of + 0] - past_u(0));
        g[of + 1] += 2 * cost_w[9] * (u[of + 1] - past_u(1));
        g[of + 2] += 2 * cost_w[10] * (u[of + 2] - past_u(2));
#if USE_EXTENDED_DYNAMICS
        g[of + 3] += 2 * cost_w[11] * u[of + 3];
#else
        g[of + 3] += 2 * cost_w[11] * (u[of + 3] - past_u(3));
#endif
        if (k != 0) {
            g[of - 4] -= 2 * cost_w[8] * (u[of + 0] - past_u(0));
            g[of - 3] -= 2 * cost_w[9] * (u[of + 1] - past_u(1));
            g[of - 2] -= 2 * cost_w[10] * (u[of + 2] - past_u(2));
#if !USE_EXTENDED_DYNAMICS
            g[of - 1] -= 2 * cost_w[11] * (u[of + 3] - past_u(3));
#endif
        }
    }

    //sum the above costs
    const double cw = cost_w[7];
    const double cbfw = cost_w[12];
    double cost = cbfw * cbf_cost + cw * c_cost + cost_w[2] * t_cost + 0.15 * yawcost + 1.0 * vcost + ucost;
    g += cbfw * cbf_cost_g + cw * c_cost_g + cost_w[2] * t_cost_g + 0.15 * yawcost_g + 1.0 * vcost_g;
    double lag_cost = 0.0;
    double contour_cost = 0.0;
    for (int k = 0; k < _n_step + 1; k++) {
        lag_cost += lag_e[k];
        contour_cost += contour_e[k];
        if (k == _n_step || k == 0) {
            cost += 2 * /*exp(-cbf_cost * 1e5) * */contour_e[k] + 0.5 * contour_e[k];
            g += 2 * (/*exp(-cbf_cost * 1e5) * */contour_e_g[k]/* - 1e5 * exp(-cbf_cost * 1e5) * contour_e[k] * cbf_cost_g*/) + 0.5 * contour_e_g[k];
            cost += cost_w[1] * lag_e[k];
            g += cost_w[1] * lag_e_g[k];
        } else {
            cost += cost_w[0] * contour_e[k];
            g += cost_w[0] * contour_e_g[k];
            cost += cost_w[1] * lag_e[k];
            g += cost_w[1] * lag_e_g[k];
        }
    }

    for (int i = 0; i < grad.size(); i++) {
        grad[i] = g[i];
    }

    instance->tmp_u_ = u;

    return cost;
}

void NominalMpcc::set_w(Matrix<double, 3 + u_dim_ + 1 + u_dim_ + 1, 1> &w) {
    cost_w = w;
}

int NominalMpcc::solve(const Matrix<double, x_dim_, 1> &state,
    const SdfMap &sdf,
    const MatrixXd &ctrl_pts,
    const MatrixXd &v_ctrl_pts,
    const MatrixXd &a_ctrl_pts,
    const double ts,
    const double len,
    const Matrix<double, u_dim_, 1> last_u,
    const Vector3d disturbance_acc,
    Matrix<double, _n_step, u_dim_> &u, 
    Matrix<double, _n_step, x_dim_> &x_predict,
    Matrix<double, _n_step + 1, 1> &t_index,
    double &solve_time) {
    const double t_min = ts * 3;
    const double t_max = ctrl_pts.rows() * ts;
    vector<double> uv(u_dim_ * _n_step + _n_step + 1);
    //initial solution
    for (int k = 0; k < _n_step; k++) {
        uv[k * u_dim_ + 0] = u(k, 0);
        uv[k * u_dim_ + 1] = u(k, 1);
        uv[k * u_dim_ + 2] = u(k, 2);
        uv[k * u_dim_ + 3] = u(k, 3);
    }
    for (int k = 0; k < _n_step; k++) {
        uv[_n_step * u_dim_ + k] = t_index[k + 1] - t_index[k];
    }
    uv[_n_step * u_dim_ + _n_step] = t_index[0];

    //boundaries on optimization variables
    vector<double> lb(uv.size()), ub(uv.size());
    for (int k = 0; k < _n_step; k++) {
        lb[k * u_dim_ + 0] = -M_PI * 1.8;
        lb[k * u_dim_ + 1] = -M_PI * 1.8;
        lb[k * u_dim_ + 2] = -M_PI * 1.8;
#if USE_EXTENDED_DYNAMICS
        lb[k * u_dim_ + 3] = -10.0;
#else
        lb[k * u_dim_ + 3] = 0.0;
#endif
        ub[k * u_dim_ + 0] = M_PI * 1.8;
        ub[k * u_dim_ + 1] = M_PI * 1.8;
        ub[k * u_dim_ + 2] = M_PI * 1.8;
#if USE_EXTENDED_DYNAMICS
        ub[k * u_dim_ + 3] = 10.0;
#else
        ub[k * u_dim_ + 3] = 0.8;
#endif
    }
    for (int k = 0; k < _n_step; k++) {
        lb[_n_step * u_dim_ + k] = 0.002 * 0.1 * (t_max - t_min) / len;
        ub[_n_step * u_dim_ + k] = 40.0 * 0.1 * (t_max - t_min) / len;
    }
    lb[_n_step * u_dim_ + _n_step] = t_index[0] + 1e-3 > t_max ? t_max : (t_index[0] + 1e-3);
    ub[_n_step * u_dim_ + _n_step] = t_max;
    opt_.set_lower_bounds(lb);
    opt_.set_upper_bounds(ub);

    for (int i = 0; i < uv.size(); i++){
        if (uv[i] > ub[i]) {
            uv[i] = ub[i];
        } else if (uv[i] < lb[i]) {
            uv[i] = lb[i];
        }
    }
    
    init_state_ = state;
    sdf_ = &sdf;
    ctrl_pts_ = ctrl_pts;
    v_ctrl_pts_ = v_ctrl_pts;
    a_ctrl_pts_ = a_ctrl_pts;
    ts_ = ts;
    t_max_ = t_max;
    last_u_ = last_u;
    disturbance_acc_ = disturbance_acc;
    for (int k = 0; k < _n_step; k++) {
        ref_u_.row(k) << 0, 0, 0, this->hover_ratio_;
    }

    double minf;
    auto beforeTime = chrono::steady_clock::now();
    try{
        flag = false;
        cnt = 0;
        //solve optimization problem
        opt_.optimize(uv, minf);
        auto afterTime = chrono::steady_clock::now();
        solve_time = chrono::duration<double>(afterTime - beforeTime).count();

        //get control inputs from optimization variables
        for (int k = 0; k < _n_step; k++) {
            u(k, 0) = uv[k * u_dim_ + 0];
            u(k, 1) = uv[k * u_dim_ + 1];
            u(k, 2) = uv[k * u_dim_ + 2];
            u(k, 3) = uv[k * u_dim_ + 3];
        }
        for (int k = 0; k < _n_step; k++) {
            x_predict.row(k) = state_[k];
        }
        t_index[0] = uv[u_dim_ * _n_step + _n_step];
        for (int k = 0; k < _n_step; k++) {
            t_index[k + 1] = t_index[k] + uv[u_dim_ * _n_step + k];
        }

        return EXIT_SUCCESS;
    } catch(exception &e) {
        auto afterTime = chrono::steady_clock::now();
        solve_time = chrono::duration<double>(afterTime - beforeTime).count();
        cerr << "nlopt failed: " << e.what() << endl;

        for (int k = 0; k < _n_step; k++) {
            u(k, 0) = uv[k * u_dim_ + 0];
            u(k, 1) = uv[k * u_dim_ + 1];
            u(k, 2) = uv[k * u_dim_ + 2];
            u(k, 3) = uv[k * u_dim_ + 3];
        }
        for (int k = 0; k < _n_step; k++) {
            x_predict.row(k) = state_[k];
        }
        t_index[0] = uv[u_dim_ * _n_step + _n_step];
        for (int k = 0; k < _n_step; k++) {
            t_index[k + 1] = t_index[k] + uv[u_dim_ * _n_step + k];
        }

        uv = tmp_u_;
        vector<double> grad(uv.size());
        VectorXd grad2(uv.size());
        flag = false;
        double cost = cost_func(uv, grad, this);
        double vub[_n_step * 3];
        double vub_g[_n_step * 3 * (_n_step * u_dim_ + _n_step + 1)];
        double vub_g2[_n_step * 3 * (_n_step * u_dim_ + _n_step + 1)];
        flag = false;
        for (int i = 0; i < uv.size(); i++) {
            double delta = 1e-8;
            vector<double> uv_new = uv;
            uv_new[i] += delta;
            vector<double> grad(uv.size());
            double cost_new = cost_func(uv_new, grad, this);
            double vub_new[_n_step * 3];
            double vub_g[_n_step * 3 * (_n_step * u_dim_ + _n_step + 1)];
            for (int j = 0; j < _n_step * 3; j++) {
                vub_g2[j * uv.size() + i] = (vub_new[j] - vub[j]) / delta;
            }
            grad2[i] = (cost_new - cost) / delta;
        }
        
        cout << "grad: " << fixed << setprecision(8) << endl << vector2Vector(grad).transpose() << endl;
        cout << "grad2: " << fixed << setprecision(8) << endl << grad2.transpose() << endl;

        return EXIT_FAILURE;
    }
}
