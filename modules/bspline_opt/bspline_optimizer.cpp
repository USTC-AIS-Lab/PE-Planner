#include <iostream>
#include <iomanip>

#include "bspline_optimizer.hpp"
#include "lbfgs/lbfgs.hpp"
#include "bspline/uniform_bspline.hpp"

using namespace std;
using namespace Eigen;

inline double BsplineOptimizer::cost_function(BsplineOptimizer *instance, const VectorXd &x, VectorXd &g) {
    for (int i = 0; i < g.rows(); i++) {
        g(i) = 0;
    }
    for (int i = 0; i < instance->ctrl_pts_->rows(); i++) {
        for (int j = 0; j < instance->ctrl_pts_->cols(); j++) {
            instance->ctrl_pts_->operator()(i, j) = x(i * 3 + j);
        }
    }
    int N = x.rows() / 3;
    double scost = 0.0;
    for (int i = 1; i < N - 1; i++) {
        double dx = x((i + 1) * 3 + 0) + x((i - 1) * 3 + 0) - 2 * x((i) * 3 + 0);
        double dy = x((i + 1) * 3 + 1) + x((i - 1) * 3 + 1) - 2 * x((i) * 3 + 1);
        double dz = x((i + 1) * 3 + 2) + x((i - 1) * 3 + 2) - 2 * x((i) * 3 + 2);
        scost += (dx * dx + dy * dy + dz * dz) * instance->lamda_s_;
        g((i + 1) * 3 + 0) += 2 * dx * instance->lamda_s_;
        g((i - 1) * 3 + 0) += 2 * dx * instance->lamda_s_;
        g((i) * 3 + 0) += 2 * dx * (-2) * instance->lamda_s_;
        g((i + 1) * 3 + 1) += 2 * dy * instance->lamda_s_;
        g((i - 1) * 3 + 1) += 2 * dy * instance->lamda_s_;
        g((i) * 3 + 1) += 2 * dy * (-2) * instance->lamda_s_;
        g((i + 1) * 3 + 2) += 2 * dz * instance->lamda_s_;
        g((i - 1) * 3 + 2) += 2 * dz * instance->lamda_s_;
        g((i) * 3 + 2) += 2 * dz * (-2) * instance->lamda_s_;
    }
    double ccost = 0.0;
    for (int i = 3; i < N - 3; i++) {
        auto sdf = instance->sdf_map_->get_dist_with_grad_trilinear(Vector3d(x(i * 3 + 0), x(i * 3 + 1), x(i * 3 + 2)));
        double dis = sdf.first;
        Vector3d grad(sdf.second(0), sdf.second(1), sdf.second(2));
        if (dis > instance->risk_dis_) {

        } else {
            ccost += pow((dis - instance->risk_dis_), 2)
                * instance->lamda_c_;
            g(i * 3) += 2 * (dis - instance->risk_dis_) * grad.x()
                * instance->lamda_c_;
            g(i * 3 + 1) += 2 * (dis - instance->risk_dis_) * grad.y()
                * instance->lamda_c_;
            g(i * 3 + 2) += 2 * (dis - instance->risk_dis_) * grad.z()
                * instance->lamda_c_;
        }
    }
    if (instance->dynobs_) {
        for (double t = 3 * instance->ts_;
            t < ((3 * instance->ts_ + 4.0) < (N * instance->ts_) ? (3 * instance->ts_ + 4.0) : (N * instance->ts_));
            t += 0.05) {
            Vector3d tmp;
            Vector4d grad2;
            int idx;
            auto p = UniformBspline::getBsplineValueFast(instance->ts_, *(instance->ctrl_pts_), t, 3, &tmp, &grad2, &idx);
            for (auto &o : *instance->dynobs_) {
                double dis;
                Vector3d grad;
                o.get_dis_ellipsoid2(p, t - 3 * instance->ts_, &dis, &grad);
                if (dis > 0.4) {

                } else {
                    // cout << dis << endl;
                    ccost += (0.4 - dis)
                        * instance->lamda_c_;
                    for (int i = 0; i < 4; i++) {
                        g((idx + i) * 3) += -grad.x() * grad2(i)
                            * instance->lamda_c_;
                        g((idx + i) * 3 + 1) += -grad.y() * grad2(i)
                            * instance->lamda_c_;
                        g((idx + i) * 3 + 2) += -grad.z() * grad2(i)
                            * instance->lamda_c_;
                    }
                }
            }
            ccost += pow(p.z() - 1.1, 2) * instance->lamda_c_ * 0.01;
            for (int i = 0; i < 4; i++) {
                g((idx + i) * 3 + 2) += 2 * (p.z() - 1.1) * grad2(i) * instance->lamda_c_ * 0.01;
            }
        }
    }
    double vcost = 0.0;
    for (int i = 0; i < N - 1; i++) {
        double vx = (x((i + 1) * 3 + 0) - x(i * 3 + 0)) / instance->ts_;
        double vy = (x((i + 1) * 3 + 1) - x(i * 3 + 1)) / instance->ts_;
        double vz = (x((i + 1) * 3 + 2) - x(i * 3 + 2)) / instance->ts_;
        if (vx * vx > pow(instance->vmax_, 2)) {
            vcost += pow(vx * vx - pow(instance->vmax_, 2), 2)
                * instance->lamda_v_;
            g(i * 3 + 0) += 2 * (vx * vx - pow(instance->vmax_, 2)) * 2 * vx
                * (-1.) / instance->ts_ * instance->lamda_v_;
            g((i + 1) * 3 + 0) += 2 * (vx * vx - pow(instance->vmax_, 2)) * 2 * vx
                / instance->ts_ * instance->lamda_v_;
        }

        if (vy * vy > pow(instance->vmax_, 2)) {
            vcost += pow((vy * vy - pow(instance->vmax_, 2)), 2)
                 * instance->lamda_v_;
            g(i * 3 + 1) += 2 * (vy * vy - pow(instance->vmax_, 2)) * 2 * vy
                * (-1.) / instance->ts_ * instance->lamda_v_;
            g((i + 1) * 3 + 1) += 2 * (vy * vy - pow(instance->vmax_, 2)) * 2 * vy
                / instance->ts_ * instance->lamda_v_;
        }

        if (vz * vz > pow(instance->vmax_, 2)) {
            vcost += pow((vz * vz - pow(instance->vmax_, 2)), 2)
                 * instance->lamda_v_;
            g(i * 3 + 2) += 2 * (vz * vz - pow(instance->vmax_, 2)) * 2 * vz
                * (-1.) / instance->ts_ * instance->lamda_v_;
            g((i + 1) * 3 + 2) += 2 * (vz * vz - pow(instance->vmax_, 2)) * 2 * vz
                / instance->ts_ * instance->lamda_v_;
        }
    }
    double acost = 0.0;
    for (int i = 0; i < N - 2; i++) {
        double ax = (x((i + 2) * 3 + 0) + x(i * 3 + 0) - 2 * x((i + 1) * 3 + 0)) 
            / pow(instance->ts_, 2);
        double ay = (x((i + 2) * 3 + 1) + x(i * 3 + 1) - 2 * x((i + 1) * 3 + 1)) 
            / pow(instance->ts_, 2);
        double az = (x((i + 2) * 3 + 2) + x(i * 3 + 2) - 2 * x((i + 1) * 3 + 2)) 
            / pow(instance->ts_, 2);
        if (ax * ax > pow(instance->amax_, 2)) {
            acost += pow((ax * ax - pow(instance->amax_, 2)), 2)
                 * instance->lamda_a_;
            g(i * 3 + 0) += 2 * (ax * ax - pow(instance->amax_, 2)) * 2 * ax
                / pow(instance->ts_, 2) * instance->lamda_a_;
            g((i + 1) * 3 + 0) += 2 * (ax * ax - pow(instance->amax_, 2)) * 2 * ax
                * (-2.) / pow(instance->ts_, 2) * instance->lamda_a_;
            g((i + 2) * 3 + 0) += 2 * (ax * ax - pow(instance->amax_, 2)) * 2 * ax
                / pow(instance->ts_, 2) * instance->lamda_a_;
        }

        if (ay * ay > pow(instance->amax_, 2)) {
            acost += pow((ay * ay - pow(instance->amax_, 2)), 2)
                 * instance->lamda_a_;
            g(i * 3 + 1) += 2 * (ay * ay - pow(instance->amax_, 2)) * 2 * ay
                / pow(instance->ts_, 2) * instance->lamda_a_;
            g((i + 1) * 3 + 1) += 2 * (ay * ay - pow(instance->amax_, 2)) * 2 * ay
                * (-2.) / pow(instance->ts_, 2) * instance->lamda_a_;
            g((i + 2) * 3 + 1) += 2 * (ay * ay - pow(instance->amax_, 2)) * 2 * ay
                / pow(instance->ts_, 2) * instance->lamda_a_;
        }

        if (az * az > pow(instance->amax_, 2)) {
            acost += pow((az * az - pow(instance->amax_, 2)), 2)
                 * instance->lamda_a_;
            g(i * 3 + 2) += 2 * (az * az - pow(instance->amax_, 2)) * 2 * az
                / pow(instance->ts_, 2) * instance->lamda_a_;
            g((i + 1) * 3 + 2) += 2 * (az * az - pow(instance->amax_, 2)) * 2 * az
                * (-2.) / pow(instance->ts_, 2) * instance->lamda_a_;
            g((i + 2) * 3 + 2) += 2 * (az * az - pow(instance->amax_, 2)) * 2 * az
                / pow(instance->ts_, 2) * instance->lamda_a_;
        }
    }
    double mean_pdis = 0.0;
    VectorXd mean_pdis_g(g.size());
    mean_pdis_g.setZero();
    // for (int i = 0; i < N - 1; i++) {
    //     double dx = x((i + 1) * 3 + 0) - x(i * 3 + 0);
    //     double dy = x((i + 1) * 3 + 1) - x(i * 3 + 1);
    //     double dz = x((i + 1) * 3 + 2) - x(i * 3 + 2);
    //     mean_pdis += (pow(dx, 2) + pow(dy, 2) + pow(dz, 2)) / (N - 1);
    //     mean_pdis_g[(i) * 3 + 0] -= 2 * dx / (N - 1);
    //     mean_pdis_g[(i) * 3 + 1] -= 2 * dy / (N - 1);
    //     mean_pdis_g[(i) * 3 + 2] -= 2 * dz / (N - 1);
    //     mean_pdis_g[(i + 1) * 3 + 0] += 2 * dx / (N - 1);
    //     mean_pdis_g[(i + 1) * 3 + 1] += 2 * dy / (N - 1);
    //     mean_pdis_g[(i + 1) * 3 + 2] += 2 * dz / (N - 1);
    // }
    // double lcost = 0.0;
    // for (int i = 0; i < N - 1; i++) {
    //     double dx = x((i + 1) * 3 + 0) - x(i * 3 + 0);
    //     double dy = x((i + 1) * 3 + 1) - x(i * 3 + 1);
    //     double dz = x((i + 1) * 3 + 2) - x(i * 3 + 2);
    //     double _c = pow(dx, 2) + pow(dy, 2) + pow(dz, 2) - mean_pdis;
    //     double c = pow(_c, 2);
    //     lcost += c * instance->lamda_l_;
    //     g((i) * 3 + 0) += 2 * _c * (-2 * dx - mean_pdis_g[(i) * 3 + 0]) * instance->lamda_l_;
    //     g((i) * 3 + 1) += 2 * _c * (-2 * dy - mean_pdis_g[(i) * 3 + 1]) * instance->lamda_l_;
    //     g((i) * 3 + 2) += 2 * _c * (-2 * dz - mean_pdis_g[(i) * 3 + 2]) * instance->lamda_l_;
    //     g((i + 1) * 3 + 0) += 2 * _c * (2 * dx - mean_pdis_g[(i + 1) * 3 + 0]) * instance->lamda_l_;
    //     g((i + 1) * 3 + 1) += 2 * _c * (2 * dy - mean_pdis_g[(i + 1) * 3 + 1]) * instance->lamda_l_;
    //     g((i + 1) * 3 + 2) += 2 * _c * (2 * dz - mean_pdis_g[(i + 1) * 3 + 2]) * instance->lamda_l_;
    // }
    for (int i = 0; i < N - 2; i++) {
        double dx = x((i + 2) * 3 + 0) - x(i * 3 + 0);
        double dy = x((i + 2) * 3 + 1) - x(i * 3 + 1);
        double dz = x((i + 2) * 3 + 2) - x(i * 3 + 2);
        mean_pdis += (pow(dx, 2) + pow(dy, 2) + pow(dz, 2)) / (N - 2);
        mean_pdis_g[(i) * 3 + 0] -= 2 * dx / (N - 2);
        mean_pdis_g[(i) * 3 + 1] -= 2 * dy / (N - 2);
        mean_pdis_g[(i) * 3 + 2] -= 2 * dz / (N - 2);
        mean_pdis_g[(i + 2) * 3 + 0] += 2 * dx / (N - 2);
        mean_pdis_g[(i + 2) * 3 + 1] += 2 * dy / (N - 2);
        mean_pdis_g[(i + 2) * 3 + 2] += 2 * dz / (N - 2);
    }
    double dlcost = 0.0;
    for (int i = 0; i < N - 2; i++) {
        double dx = x((i + 2) * 3 + 0) - x(i * 3 + 0);
        double dy = x((i + 2) * 3 + 1) - x(i * 3 + 1);
        double dz = x((i + 2) * 3 + 2) - x(i * 3 + 2);
        double _c = pow(dx, 2) + pow(dy, 2) + pow(dz, 2) - mean_pdis;
        double c = pow(_c, 2);
        dlcost += c * instance->lamda_dl_;
        g((i) * 3 + 0) += 2 * _c * (-2 * dx - mean_pdis_g[(i) * 3 + 0]) * instance->lamda_dl_;
        g((i) * 3 + 1) += 2 * _c * (-2 * dy - mean_pdis_g[(i) * 3 + 1]) * instance->lamda_dl_;
        g((i) * 3 + 2) += 2 * _c * (-2 * dz - mean_pdis_g[(i) * 3 + 2]) * instance->lamda_dl_;
        g((i + 2) * 3 + 0) += 2 * _c * (2 * dx - mean_pdis_g[(i + 2) * 3 + 0]) * instance->lamda_dl_;
        g((i + 2) * 3 + 1) += 2 * _c * (2 * dy - mean_pdis_g[(i + 2) * 3 + 1]) * instance->lamda_dl_;
        g((i + 2) * 3 + 2) += 2 * _c * (2 * dz - mean_pdis_g[(i + 2) * 3 + 2]) * instance->lamda_dl_;
    }
    double lcost = instance->lamda_l_ * mean_pdis * (N - 2);
    g += instance->lamda_l_ * mean_pdis_g * (N - 2);
    double ecost = 0.0;
    {
        {
            double dx = (x(0) + 4 * x(3 + 0) + x(2 * 3 + 0)) / 6.0 - instance->start_p_(0);
            double dy = (x(1) + 4 * x(3 + 1) + x(2 * 3 + 1)) / 6.0 - instance->start_p_(1);
            double dz = (x(2) + 4 * x(3 + 2) + x(2 * 3 + 2)) / 6.0 - instance->start_p_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ep_;
            g(0) += 2 * dx / 6.0 * instance->lamda_ep_;
            g(1) += 2 * dy / 6.0 * instance->lamda_ep_;
            g(2) += 2 * dz / 6.0 * instance->lamda_ep_;
            g(3) += 8 * dx / 6.0 * instance->lamda_ep_;
            g(4) += 8 * dy / 6.0 * instance->lamda_ep_;
            g(5) += 8 * dz / 6.0 * instance->lamda_ep_;
            g(6) += 2 * dx / 6.0 * instance->lamda_ep_;
            g(7) += 2 * dy / 6.0 * instance->lamda_ep_;
            g(8) += 2 * dz / 6.0 * instance->lamda_ep_;
        }
        {
            double dx = (x((N - 3) * 3 + 0) + 4 * x((N - 2) * 3 + 0) + x((N - 1) * 3 + 0)) / 6.0 - instance->end_p_(0);
            double dy = (x((N - 3) * 3 + 1) + 4 * x((N - 2) * 3 + 1) + x((N - 1) * 3 + 1)) / 6.0 - instance->end_p_(1);
            double dz = (x((N - 3) * 3 + 2) + 4 * x((N - 2) * 3 + 2) + x((N - 1) * 3 + 2)) / 6.0 - instance->end_p_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ep_;
            g((N - 3) * 3 + 0) += 2 * dx / 6.0 * instance->lamda_ep_;
            g((N - 3) * 3 + 1) += 2 * dy / 6.0 * instance->lamda_ep_;
            g((N - 3) * 3 + 2) += 2 * dz / 6.0 * instance->lamda_ep_;
            g((N - 2) * 3 + 0) += 8 * dx / 6.0 * instance->lamda_ep_;
            g((N - 2) * 3 + 1) += 8 * dy / 6.0 * instance->lamda_ep_;
            g((N - 2) * 3 + 2) += 8 * dz / 6.0 * instance->lamda_ep_;
            g((N - 1) * 3 + 0) += 2 * dx / 6.0 * instance->lamda_ep_;
            g((N - 1) * 3 + 1) += 2 * dy / 6.0 * instance->lamda_ep_;
            g((N - 1) * 3 + 2) += 2 * dz / 6.0 * instance->lamda_ep_;
        }
        {
            double dx = (-x(0) + x(2 * 3 + 0)) / 2.0 / instance->ts_ - instance->start_v_(0);
            double dy = (-x(1) + x(2 * 3 + 1)) / 2.0 / instance->ts_ - instance->start_v_(1);
            double dz = (-x(2) + x(2 * 3 + 2)) / 2.0 / instance->ts_ - instance->start_v_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ev_;
            g(0) -= 2 * dx / 2.0 / instance->ts_ * instance->lamda_ev_;
            g(1) -= 2 * dy / 2.0 / instance->ts_ * instance->lamda_ev_;
            g(2) -= 2 * dz / 2.0 / instance->ts_ * instance->lamda_ev_;
            g(6) += 2 * dx / 2.0 / instance->ts_ * instance->lamda_ev_;
            g(7) += 2 * dy / 2.0 / instance->ts_ * instance->lamda_ev_;
            g(8) += 2 * dz / 2.0 / instance->ts_ * instance->lamda_ev_;
        }
        {
            double dx = (-x((N - 3) * 3 + 0) + x((N - 1) * 3 + 0)) / 2.0 / instance->ts_ - instance->end_v_(0);
            double dy = (-x((N - 3) * 3 + 1) + x((N - 1) * 3 + 1)) / 2.0 / instance->ts_ - instance->end_v_(1);
            double dz = (-x((N - 3) * 3 + 2) + x((N - 1) * 3 + 2)) / 2.0 / instance->ts_ - instance->end_v_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ev_;
            g((N - 3) * 3 + 0) -= 2 * dx / 2.0 / instance->ts_ * instance->lamda_ev_;
            g((N - 3) * 3 + 1) -= 2 * dy / 2.0 / instance->ts_ * instance->lamda_ev_;
            g((N - 3) * 3 + 2) -= 2 * dz / 2.0 / instance->ts_ * instance->lamda_ev_;
            g((N - 1) * 3 + 0) += 2 * dx / 2.0 / instance->ts_ * instance->lamda_ev_;
            g((N - 1) * 3 + 1) += 2 * dy / 2.0 / instance->ts_ * instance->lamda_ev_;
            g((N - 1) * 3 + 2) += 2 * dz / 2.0 / instance->ts_ * instance->lamda_ev_;
        }
        {
            double t2 = instance->ts_ * instance->ts_;
            double dx = (x(0) - 2 * x(3 + 0) + x(2 * 3 + 0)) / t2 - instance->start_a_(0);
            double dy = (x(1) - 2 * x(3 + 1) + x(2 * 3 + 1)) / t2 - instance->start_a_(1);
            double dz = (x(2) - 2 * x(3 + 2) + x(2 * 3 + 2)) / t2 - instance->start_a_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ea_;
            g(0) += 2 * dx / t2 * instance->lamda_ea_;
            g(1) += 2 * dy / t2 * instance->lamda_ea_;
            g(2) += 2 * dz / t2 * instance->lamda_ea_;
            g(3) -= 4 * dx / t2 * instance->lamda_ea_;
            g(4) -= 4 * dy / t2 * instance->lamda_ea_;
            g(5) -= 4 * dz / t2 * instance->lamda_ea_;
            g(6) += 2 * dx / t2 * instance->lamda_ea_;
            g(7) += 2 * dy / t2 * instance->lamda_ea_;
            g(8) += 2 * dz / t2 * instance->lamda_ea_;
        }
        {
            double t2 = instance->ts_ * instance->ts_;
            double dx = (x((N - 3) * 3 + 0) - 2 * x((N - 2) * 3 + 0) + x((N - 1) * 3 + 0)) / t2 - instance->end_a_(0);
            double dy = (x((N - 3) * 3 + 1) - 2 * x((N - 2) * 3 + 1) + x((N - 1) * 3 + 1)) / t2 - instance->end_a_(1);
            double dz = (x((N - 3) * 3 + 2) - 2 * x((N - 2) * 3 + 2) + x((N - 1) * 3 + 2)) / t2 - instance->end_a_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ea_;
            g((N - 3) * 3 + 0) += 2 * dx / t2 * instance->lamda_ea_;
            g((N - 3) * 3 + 1) += 2 * dy / t2 * instance->lamda_ea_;
            g((N - 3) * 3 + 2) += 2 * dz / t2 * instance->lamda_ea_;
            g((N - 2) * 3 + 0) += -4 * dx / t2 * instance->lamda_ea_;
            g((N - 2) * 3 + 1) += -4 * dy / t2 * instance->lamda_ea_;
            g((N - 2) * 3 + 2) += -4 * dz / t2 * instance->lamda_ea_;
            g((N - 1) * 3 + 0) += 2 * dx / t2 * instance->lamda_ea_;
            g((N - 1) * 3 + 1) += 2 * dy / t2 * instance->lamda_ea_;
            g((N - 1) * 3 + 2) += 2 * dz / t2 * instance->lamda_ea_;
        }
    }

    instance->cnt_++;
    return scost + ccost + vcost + acost + lcost + dlcost + ecost;
}

int BsplineOptimizer::optimize(MatrixXd &ctrl_pts
    , const Vector3d &start_p, const Vector3d &start_v, const Vector3d &start_a
    , const Vector3d &end_p, const Vector3d &end_v, const Vector3d &end_a
    , const SdfMap *const sdf, const double ts
    , const double lamda_s, const double lamda_c, const double lamda_v
    , const double lamda_a, const double lamda_l, const double lamda_dl
    , const double lamda_ep, const double lamda_ev, const double lamda_ea
    , const double risk_dis, const double vmax
    , const double amax, vector<tuple<double, double, double, double, double, double>> &cost_history
    , const vector<DynObs> *dynobs) {
    VectorXd x = VectorXd::Zero(ctrl_pts.rows() * ctrl_pts.cols());
    for (int i = 0; i < ctrl_pts.rows(); i++) {
        for (int j = 0; j < ctrl_pts.cols(); j++) {
            x(i * 3 + j) = ctrl_pts(i, j);
        }
    }
    sdf_map_ = sdf;
    dynobs_ = dynobs;
    ts_ = ts;
    lamda_s_ = lamda_s;
    lamda_c_ = lamda_c;
    lamda_v_ = lamda_v;
    lamda_a_ = lamda_a;
    lamda_l_ = lamda_l;
    lamda_dl_ = lamda_dl;
    lamda_ep_ = lamda_ep;
    lamda_ev_ = lamda_ev;
    lamda_ea_ = lamda_ea;
    risk_dis_ = risk_dis;
    start_p_ = start_p;
    start_v_ = start_v;
    start_a_ = start_a;
    end_p_ = end_p;
    end_v_ = end_v;
    end_a_ = end_a;
    vmax_ = vmax;
    amax_ = amax;
    ctrl_pts_ = &ctrl_pts;

    lbfgs::lbfgs_parameter_t params;
    params.g_epsilon = 1.0e-10;
    params.past = 3;
    params.delta = 1.0e-6;
    params.mem_size = 8;
    params.max_linesearch = 128;

    double final_cost = 0.0;
    cost_history_.clear();
    cnt_ = 0;
    time_ = chrono::steady_clock::now();
    int ret = lbfgs::lbfgs_optimize(x,
        final_cost,
        (lbfgs::lbfgs_evaluate_t)cost_function,
        nullptr,
        nullptr,
        this,
        params);
    double spend = chrono::duration<double>(chrono::steady_clock::now() - time_).count();
    cost_history = cost_history_;

    // cout << fixed << setprecision(4)
    //     << "L-BFGS Optimization Returned: " << ret << endl
    //     << "B-spline optimization spend " << spend * 1e3 << " ms" << endl
    //     << "Minimized Cost: " << final_cost << endl;
    
    for (int i = 0; i < ctrl_pts.rows(); i++) {
        for (int j = 0; j < ctrl_pts.cols(); j++) {
            ctrl_pts(i, j) = x(i * 3 + j);
        }
    }

    return ret;
}
