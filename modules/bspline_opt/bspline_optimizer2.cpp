#include "bspline_optimizer2.hpp"

#include <unistd.h>


BsplineOptimizer2::BsplineOptimizer2() {

}

BsplineOptimizer2::~BsplineOptimizer2() {

}

int BsplineOptimizer2::optimize(MatrixXd &ctrl_pts
        , const Vector3d &start_p, const Vector3d &start_v
        , const Vector3d &end_p, const Vector3d &end_v
        , const SdfMap *const sdf, const double ts
        , const double lamda_s, const double lamda_c, const double lamda_v
        , const double lamda_a, const double lamda_l, const double lamda_ep, const double lamda_ev
        , const double risk_dis, const double vmax
        , const double amax, vector<tuple<double, double, double, double, double, double>> &cost_history) {
    nlopt::opt opt("LD_LBFGS", ctrl_pts.size());
    opt.set_min_objective((nlopt::vfunc)&BsplineOptimizer2::cost_func, this);
    opt.set_xtol_rel(1e-12);
    opt.set_ftol_rel(1e-9);
    opt.set_maxeval(1000);
    opt.set_vector_storage(8);

    // vector<double> lb, ub;
    // for (int i = 0; i < ctrl_pts.size(); i++) {
    //     lb.push_back(-1e6);
    //     ub.push_back(1e6);
    // }
    // opt.set_lower_bounds(lb);
    // opt.set_upper_bounds(ub);

    sdf_map_ = sdf;
    ts_ = ts;
    lamda_s_ = lamda_s;
    lamda_c_ = lamda_c;
    lamda_v_ = lamda_v;
    lamda_a_ = lamda_a;
    lamda_l_ = lamda_l;
    lamda_ep_ = lamda_ep;
    lamda_ev_ = lamda_ev;
    risk_dis_ = risk_dis;
    start_p_ = start_p;
    start_v_ = start_v;
    end_p_ = end_p;
    end_v_ = end_v;
    vmax_ = vmax;
    amax_ = amax;

    cost_history_.clear();

    vector<double> uv(ctrl_pts.size());
    for (int i = 0; i < ctrl_pts.rows(); i++) {
        for (int j = 0; j < ctrl_pts.cols(); j++) {
            uv[i * 3 + j] = ctrl_pts(i, j);
        }
    }

    double minf = 0.0;

    try{
        cnt_ = 0;
        time_ = chrono::steady_clock::now();
        opt.optimize(uv, minf);
        auto afterTime = chrono::steady_clock::now();
        double solve_time = chrono::duration<double>(afterTime - time_).count();
        for (int i = 0; i < ctrl_pts.rows(); i++) {
            for (int j = 0; j < ctrl_pts.cols(); j++) {
                ctrl_pts(i, j) = uv[i * 3 + j];
            }
        }
        cost_history = cost_history_;
        cout << "B-spline optimization spend " << solve_time * 1e3 << " ms" << endl;
        cout << "Minimized Cost: " << minf << endl;
        return 0;
    } catch(exception &e) {
        cerr << "nlopt failed: " << e.what() << endl;
        return 1;
    }
}

inline double BsplineOptimizer2::cost_func(const vector<double> &u, vector<double> &g, BsplineOptimizer2 *instance) {
    for (int i = 0; i < g.size(); i++) {
        g[i] = 0;
    }
    int N = u.size() / 3;
    double scost = 0.0;
    for (int i = 1; i < N - 1; i++) {
        double dx = u[(i + 1) * 3 + 0] + u[(i - 1) * 3 + 0] - 2 * u[(i) * 3 + 0];
        double dy = u[(i + 1) * 3 + 1] + u[(i - 1) * 3 + 1] - 2 * u[(i) * 3 + 1];
        double dz = u[(i + 1) * 3 + 2] + u[(i - 1) * 3 + 2] - 2 * u[(i) * 3 + 2];
        scost += (dx * dx + dy * dy + dz * dz) * instance->lamda_s_;
        g[(i + 1) * 3 + 0] += 2 * dx * instance->lamda_s_;
        g[(i - 1) * 3 + 0] += 2 * dx * instance->lamda_s_;
        g[(i) * 3 + 0] += 2 * dx * (-2) * instance->lamda_s_;
        g[(i + 1) * 3 + 1] += 2 * dy * instance->lamda_s_;
        g[(i - 1) * 3 + 1] += 2 * dy * instance->lamda_s_;
        g[(i) * 3 + 1] += 2 * dy * (-2) * instance->lamda_s_;
        g[(i + 1) * 3 + 2] += 2 * dz * instance->lamda_s_;
        g[(i - 1) * 3 + 2] += 2 * dz * instance->lamda_s_;
        g[(i) * 3 + 2] += 2 * dz * (-2) * instance->lamda_s_;
    }
    double ccost = 0.0;
    for (int i = 0; i < N; i++) {
        auto sdf = instance->sdf_map_->get_dist_with_grad_trilinear(Vector3d(u[i * 3 + 0], u[i * 3 + 1], u[i * 3 + 2]));
        double dis = sdf.first;
        Vector3d grad(sdf.second(0), sdf.second(1), sdf.second(2));
        if (dis > instance->risk_dis_) {

        } else {
            ccost += pow((dis - instance->risk_dis_), 2)
                * instance->lamda_c_;
            g[i * 3] += 2 * (dis - instance->risk_dis_) * grad.x()
                * instance->lamda_c_;
            g[i * 3 + 1] += 2 * (dis - instance->risk_dis_) * grad.y()
                * instance->lamda_c_;
            g[i * 3 + 2] += 2 * (dis - instance->risk_dis_) * grad.z()
                * instance->lamda_c_;
        }
    }
    double vcost = 0.0;
    for (int i = 0; i < N - 1; i++) {
        double vx = (u[(i + 1) * 3 + 0] - u[i * 3 + 0]) / instance->ts_;
        double vy = (u[(i + 1) * 3 + 1] - u[i * 3 + 1]) / instance->ts_;
        double vz = (u[(i + 1) * 3 + 2] - u[i * 3 + 2]) / instance->ts_;
        if (vx * vx > pow(instance->vmax_, 2)) {
            vcost += pow(vx * vx - pow(instance->vmax_, 2), 2)
                * instance->lamda_v_;
            g[i * 3 + 0] += 2 * (vx * vx - pow(instance->vmax_, 2)) * 2 * vx
                * (-1.) / instance->ts_ * instance->lamda_v_;
            g[(i + 1) * 3 + 0] += 2 * (vx * vx - pow(instance->vmax_, 2)) * 2 * vx
                / instance->ts_ * instance->lamda_v_;
        }

        if (vy * vy > pow(instance->vmax_, 2)) {
            vcost += pow((vy * vy - pow(instance->vmax_, 2)), 2)
                 * instance->lamda_v_;
            g[i * 3 + 1] += 2 * (vy * vy - pow(instance->vmax_, 2)) * 2 * vy
                * (-1.) / instance->ts_ * instance->lamda_v_;
            g[(i + 1) * 3 + 1] += 2 * (vy * vy - pow(instance->vmax_, 2)) * 2 * vy
                / instance->ts_ * instance->lamda_v_;
        }

        if (vz * vz > pow(instance->vmax_, 2)) {
            vcost += pow((vz * vz - pow(instance->vmax_, 2)), 2)
                 * instance->lamda_v_;
            g[i * 3 + 2] += 2 * (vz * vz - pow(instance->vmax_, 2)) * 2 * vz
                * (-1.) / instance->ts_ * instance->lamda_v_;
            g[(i + 1) * 3 + 2] += 2 * (vz * vz - pow(instance->vmax_, 2)) * 2 * vz
                / instance->ts_ * instance->lamda_v_;
        }
    }
    double acost = 0.0;
    for (int i = 0; i < N - 2; i++) {
        double ax = (u[(i + 2) * 3 + 0] + u[i * 3 + 0] - 2 * u[(i + 1) * 3 + 0]) 
            / pow(instance->ts_, 2);
        double ay = (u[(i + 2) * 3 + 1] + u[i * 3 + 1] - 2 * u[(i + 1) * 3 + 1]) 
            / pow(instance->ts_, 2);
        double az = (u[(i + 2) * 3 + 2] + u[i * 3 + 2] - 2 * u[(i + 1) * 3 + 2]) 
            / pow(instance->ts_, 2);
        if (ax * ax > pow(instance->amax_, 2)) {
            acost += pow((ax * ax - pow(instance->amax_, 2)), 2)
                 * instance->lamda_a_;
            g[i * 3 + 0] += 2 * (ax * ax - pow(instance->amax_, 2)) * 2 * ax
                / pow(instance->ts_, 2) * instance->lamda_a_;
            g[(i + 1) * 3 + 0] += 2 * (ax * ax - pow(instance->amax_, 2)) * 2 * ax
                * (-2.) / pow(instance->ts_, 2) * instance->lamda_a_;
            g[(i + 2) * 3 + 0] += 2 * (ax * ax - pow(instance->amax_, 2)) * 2 * ax
                / pow(instance->ts_, 2) * instance->lamda_a_;
        }

        if (ay * ay > pow(instance->amax_, 2)) {
            acost += pow((ay * ay - pow(instance->amax_, 2)), 2)
                 * instance->lamda_a_;
            g[i * 3 + 1] += 2 * (ay * ay - pow(instance->amax_, 2)) * 2 * ay
                / pow(instance->ts_, 2) * instance->lamda_a_;
            g[(i + 1) * 3 + 1] += 2 * (ay * ay - pow(instance->amax_, 2)) * 2 * ay
                * (-2.) / pow(instance->ts_, 2) * instance->lamda_a_;
            g[(i + 2) * 3 + 1] += 2 * (ay * ay - pow(instance->amax_, 2)) * 2 * ay
                / pow(instance->ts_, 2) * instance->lamda_a_;
        }

        if (az * az > pow(instance->amax_, 2)) {
            acost += pow((az * az - pow(instance->amax_, 2)), 2)
                 * instance->lamda_a_;
            g[i * 3 + 2] += 2 * (az * az - pow(instance->amax_, 2)) * 2 * az
                / pow(instance->ts_, 2) * instance->lamda_a_;
            g[(i + 1) * 3 + 2] += 2 * (az * az - pow(instance->amax_, 2)) * 2 * az
                * (-2.) / pow(instance->ts_, 2) * instance->lamda_a_;
            g[(i + 2) * 3 + 2] += 2 * (az * az - pow(instance->amax_, 2)) * 2 * az
                / pow(instance->ts_, 2) * instance->lamda_a_;
        }
    }
    double mean_pdis = 0.0; //控制点平均距离
    VectorXd mean_pdis_g(g.size());
    mean_pdis_g.setZero();
    for (int i = 0; i < N - 1; i++) {
        double dx = u[(i + 1) * 3 + 0] - u[i * 3 + 0];
        double dy = u[(i + 1) * 3 + 1] - u[i * 3 + 1];
        double dz = u[(i + 1) * 3 + 2] - u[i * 3 + 2];
        mean_pdis += (pow(dx, 2) + pow(dy, 2) + pow(dz, 2)) / (N - 1);
        mean_pdis_g[(i) * 3 + 0] -= 2 * dx / (N - 1);
        mean_pdis_g[(i) * 3 + 1] -= 2 * dy / (N - 1);
        mean_pdis_g[(i) * 3 + 2] -= 2 * dz / (N - 1);
        mean_pdis_g[(i + 1) * 3 + 0] += 2 * dx / (N - 1);
        mean_pdis_g[(i + 1) * 3 + 1] += 2 * dy / (N - 1);
        mean_pdis_g[(i + 1) * 3 + 2] += 2 * dz / (N - 1);
    }
    double lcost = 0.0;
    for (int i = 0; i < N - 1; i++) {
        double dx = u[(i + 1) * 3 + 0] - u[i * 3 + 0];
        double dy = u[(i + 1) * 3 + 1] - u[i * 3 + 1];
        double dz = u[(i + 1) * 3 + 2] - u[i * 3 + 2];
        double _c = pow(dx, 2) + pow(dy, 2) + pow(dz, 2) - mean_pdis;
        double c = pow(_c, 2);
        lcost += c * instance->lamda_l_;
        g[(i) * 3 + 0] += 2 * _c * (-2 * dx - mean_pdis_g[(i) * 3 + 0]) * instance->lamda_l_;
        g[(i) * 3 + 1] += 2 * _c * (-2 * dy - mean_pdis_g[(i) * 3 + 1]) * instance->lamda_l_;
        g[(i) * 3 + 2] += 2 * _c * (-2 * dz - mean_pdis_g[(i) * 3 + 2]) * instance->lamda_l_;
        g[(i + 1) * 3 + 0] += 2 * _c * (2 * dx - mean_pdis_g[(i + 1) * 3 + 0]) * instance->lamda_l_;
        g[(i + 1) * 3 + 1] += 2 * _c * (2 * dy - mean_pdis_g[(i + 1) * 3 + 1]) * instance->lamda_l_;
        g[(i + 1) * 3 + 2] += 2 * _c * (2 * dz - mean_pdis_g[(i + 1) * 3 + 2]) * instance->lamda_l_;
    }

    //端代价
    double ecost = 0.0;
    {
        {
            double dx = (u[0] + 4 * u[3 + 0] + u[2 * 3 + 0]) / 6.0 - instance->start_p_(0);
            double dy = (u[1] + 4 * u[3 + 1] + u[2 * 3 + 1]) / 6.0 - instance->start_p_(1);
            double dz = (u[2] + 4 * u[3 + 2] + u[2 * 3 + 2]) / 6.0 - instance->start_p_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ep_;
            g[0] += 2 * dx / 6.0 * instance->lamda_ep_;
            g[1] += 2 * dy / 6.0 * instance->lamda_ep_;
            g[2] += 2 * dz / 6.0 * instance->lamda_ep_;
            g[3] += 8 * dx / 6.0 * instance->lamda_ep_;
            g[4] += 8 * dy / 6.0 * instance->lamda_ep_;
            g[5] += 8 * dz / 6.0 * instance->lamda_ep_;
            g[6] += 2 * dx / 6.0 * instance->lamda_ep_;
            g[7] += 2 * dy / 6.0 * instance->lamda_ep_;
            g[8] += 2 * dz / 6.0 * instance->lamda_ep_;
        }
        {
            double dx = (u[(N - 3) * 3 + 0] + 4 * u[(N - 2) * 3 + 0] + u[(N - 1) * 3 + 0]) / 6.0 - instance->end_p_(0);
            double dy = (u[(N - 3) * 3 + 1] + 4 * u[(N - 2) * 3 + 1] + u[(N - 1) * 3 + 1]) / 6.0 - instance->end_p_(1);
            double dz = (u[(N - 3) * 3 + 2] + 4 * u[(N - 2) * 3 + 2] + u[(N - 1) * 3 + 2]) / 6.0 - instance->end_p_(2);
            ecost += (dx * dx + dy * dy + dz * dz) * instance->lamda_ep_;
            g[(N - 3) * 3 + 0] += 2 * dx / 6.0 * instance->lamda_ep_;
            g[(N - 3) * 3 + 1] += 2 * dy / 6.0 * instance->lamda_ep_;
            g[(N - 3) * 3 + 2] += 2 * dz / 6.0 * instance->lamda_ep_;
            g[(N - 2) * 3 + 0] += 8 * dx / 6.0 * instance->lamda_ep_;
            g[(N - 2) * 3 + 1] += 8 * dy / 6.0 * instance->lamda_ep_;
            g[(N - 2) * 3 + 2] += 8 * dz / 6.0 * instance->lamda_ep_;
            g[(N - 1) * 3 + 0] += 2 * dx / 6.0 * instance->lamda_ep_;
            g[(N - 1) * 3 + 1] += 2 * dy / 6.0 * instance->lamda_ep_;
            g[(N - 1) * 3 + 2] += 2 * dz / 6.0 * instance->lamda_ep_;
        }
    }

    instance->cnt_++;
    return scost + ccost + vcost + acost + lcost + ecost;
}