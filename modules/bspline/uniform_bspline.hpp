#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace Eigen;

class UniformBspline {
public:
    UniformBspline() {}

    static void parameter2Bspline(const double ts, const vector<Vector3d> &poss
        , const vector<Vector3d> &vels, const vector<Vector3d> &accs
        , MatrixXd &ctrl_pts);
    static void parameter2Bspline(const double ts, const vector<double> &poss
        , const vector<double> &vels
        , MatrixXd &ctrl_pts);
    static Vector3d getBsplineValue(const double ts, const MatrixXd &ctrl_pts, double t, int degree);
    static double basisFunction(const int k, const int d, const double u, const double ts);
    static MatrixXd getDerivativeCtrlPts(const MatrixXd &ctrl_pts, double ts);
    template <class T> static T getBsplineValueFast(const double ts, const MatrixXd &ctrl_pts, double t, int degree, T *grad = nullptr, Vector4d *grad2=nullptr, int *idx=nullptr);
    static double getArcLength(double ts, const MatrixXd &deriv_ctrl_pts, double dt, double start_t, double end_t, vector<double> *arc_lengths);
};

inline void UniformBspline::parameter2Bspline(const double ts, const vector<Vector3d> &poss
    , const vector<Vector3d> &vels, const vector<Vector3d> &accs
    , MatrixXd &ctrl_pts) {
    int K = poss.size();
    Vector3d prow(1, 4, 1), vrow(-1, 0, 1), arow(1, -2, 1);
    MatrixXd coefs = MatrixXd::Zero(K + 4, K + 2);
    for (int i = 0; i < K; i++) {
        coefs.block(i, i, 1, 3) = 1. / 6. * prow.transpose();
    }
    coefs.block(K, 0, 1, 3) = 1. / 2. / ts * vrow.transpose();
    coefs.block(K + 1, K - 1, 1, 3) = 1. / 2. / ts * vrow.transpose();
    coefs.block(K + 2, 0, 1, 3) = 1. / ts / ts * arow.transpose();
    coefs.block(K + 3, K - 1, 1, 3) = 1. / ts / ts * arow.transpose();
    VectorXd bx(K + 4), by(K + 4), bz(K + 4);
    for (int i = 0; i < K; i++) {
        bx(i) = poss[i](0);
        by(i) = poss[i](1);
        bz(i) = poss[i](2);
    }
    bx(K) = vels[0](0);
    by(K) = vels[0](1);
    bz(K) = vels[0](2);
    bx(K + 1) = vels[1](0);
    by(K + 1) = vels[1](1);
    bz(K + 1) = vels[1](2);
    bx(K + 2) = accs[0](0);
    by(K + 2) = accs[0](1);
    bz(K + 2) = accs[0](2);
    bx(K + 3) = accs[1](0);
    by(K + 3) = accs[1](1);
    bz(K + 3) = accs[1](2);

    auto time = chrono::steady_clock::now();
    VectorXd px = coefs.colPivHouseholderQr().solve(bx);
    VectorXd py = coefs.colPivHouseholderQr().solve(by);
    VectorXd pz = coefs.colPivHouseholderQr().solve(bz);
    ctrl_pts.resize(K + 2, 3);
    ctrl_pts.col(0) = px;
    ctrl_pts.col(1) = py;
    ctrl_pts.col(2) = pz;
    // cout << "B-spline parameterization spend " << chrono::duration<double>(chrono::steady_clock::now() - time).count() * 1e3 << " ms" << endl;
}

inline void UniformBspline::parameter2Bspline(const double ts, const vector<double> &poss
    , const vector<double> &vels
    , MatrixXd &ctrl_pts) {
    int K = poss.size();
    Vector3d prow(1, 4, 1), vrow(-1, 0, 1);
    MatrixXd coefs = MatrixXd::Zero(K + 2, K + 2);
    for (int i = 0; i < K; i++) {
        coefs.block(i, i, 1, 3) = 1. / 6. * prow.transpose();
    }
    coefs.block(K, 0, 1, 3) = 1. / 2. / ts * vrow.transpose();
    coefs.block(K + 1, K - 1, 1, 3) = 1. / 2. / ts * vrow.transpose();
    VectorXd b(K + 2);
    for (int i = 0; i < K; i++) {
        b(i) = poss[i];
    }
    b(K) = vels[0];
    b(K + 1) = vels[1];

    auto time = chrono::steady_clock::now();
    VectorXd p = coefs.colPivHouseholderQr().solve(b);
    ctrl_pts.resize(K + 2, 3);
    ctrl_pts.col(0) = p;
    cout << "B-spline parameterization spend " << chrono::duration<double>(chrono::steady_clock::now() - time).count() * 1e3 << " ms" << endl;
}

inline Vector3d UniformBspline::getBsplineValue(const double ts, const MatrixXd &ctrl_pts, double t, int degree) {
    if (t < degree * ts - 1e-6 || t > ctrl_pts.rows() * ts + 1e-6) {
        return Vector3d(INFINITY, INFINITY, INFINITY);
    }
    int time_idx = 3;
    while (t > (time_idx + 1) * ts) {
        time_idx++;
    }
    Vector3d ret = Vector3d::Zero();
    for (int i = 0; i < ctrl_pts.rows(); i++) {
        auto tmp = basisFunction(i, degree, t, ts);
        ret += ctrl_pts.row(i) * tmp;
        // cout << fixed << setprecision(32) << tmp << " ";
    }
    // cout << endl;
    return ret;
}

inline double UniformBspline::basisFunction(const int k, const int d, const double u, const double ts) {
    if (d == 0) {
        if (u > k * ts && u <= (k + 1) * ts) {
            return 1;
        } else {
            return 0;
        }
    } else {
        if (u > k * ts && u <= (k + d + 1) * ts) {
            return (u - k * ts) / (d * ts) * basisFunction(k, d - 1, u, ts)
                + ((k + d + 1) * ts - u) / (d * ts) * basisFunction(k + 1, d - 1, u, ts);
        } else {
            return 0;
        }
    }
}

inline MatrixXd UniformBspline::getDerivativeCtrlPts(const MatrixXd &ctrl_pts, double ts) {
    MatrixXd derivative_pts = MatrixXd::Zero(ctrl_pts.rows() - 1, ctrl_pts.cols());
    for (int i = 0; i < ctrl_pts.rows() - 1; i++) {
        Vector3d q1 = ctrl_pts.row(i + 1);
        Vector3d q0 = ctrl_pts.row(i);
        derivative_pts.row(i) = (q1 - q0) / ts;
    }
    return derivative_pts;
}

template <class T>
inline T UniformBspline::getBsplineValueFast(const double ts, const MatrixXd &ctrl_pts, double t, int degree, T *grad, Vector4d *grad2, int *idx) {
    if (t < degree * ts - 1e-6 || t > ctrl_pts.rows() * ts + 1e-6) {
        throw "getBsplineValueFast parameter error";
        exit(0);
        // return Vector3d(INFINITY, INFINITY, INFINITY);
    }
    if (degree < 1 || degree > 3) {
        throw "getBsplineValueFast parameter error";
        exit(0);
        // return Vector3d(INFINITY, INFINITY, INFINITY);
    }

    T ret(ctrl_pts.cols());

    int k = t / ts;
    if (k < degree) {
        k = degree;
    } else if (k > ctrl_pts.rows() - 1) {
        k = ctrl_pts.rows() - 1;
    }
    double x = (t - k * ts) / ts;
    double x_g = 1 / ts;
    double ix = ((k + 1) * ts - t) / ts;
    double ix_g = -1 / ts;
    double w0, w1, w2, w3;
    double w0_g, w1_g, w2_g, w3_g;

    switch(degree) {
    case 1:
        w0 = ix;
        w1 = x;
        ret = w0 * ctrl_pts.row(k - 1) + w1 * ctrl_pts.row(k);
        if (grad) {
            w0_g = ix_g;
            w1_g = x_g;
            *grad = w0_g * ctrl_pts.row(k - 1) + w1_g * ctrl_pts.row(k);
        }
        break;
    case 2:
        w0 = 0.5 * ix * ix;
        w1 = 0.5 + x - x * x;
        w2 = 0.5 * x * x;
        ret = w0 * ctrl_pts.row(k - 2) + w1 * ctrl_pts.row(k - 1) + w2 * ctrl_pts.row(k);
        if (grad) {
            w0_g = ix * ix_g;
            w1_g = x_g - 2 * x * x_g;
            w2_g = x * x_g;
            *grad = w0_g * ctrl_pts.row(k - 2) + w1_g * ctrl_pts.row(k - 1) + w2_g * ctrl_pts.row(k);
        }
        break;
    case 3:
        w0 = ix * ix * ix / 6.0;
        w1 = 1 / 6.0 + 0.5 * (ix + ix * ix - ix * ix * ix);
        w2 = 1 / 6.0 + 0.5 * (x + x * x - x * x * x);
        w3 = x * x * x / 6.0;
        // cout << fixed << setprecision(32) << w0 << " " << w1 << " " << w2 << " " << w3 << endl;
        ret = w0 * ctrl_pts.row(k - 3) + w1 * ctrl_pts.row(k - 2) + w2 * ctrl_pts.row(k - 1) + w3 * ctrl_pts.row(k);
        if (grad) {
            w0_g = 0.5 * ix * ix * ix_g;
            w1_g = 0.5 * (ix_g + 2 * ix * ix_g - 3 * ix * ix * ix_g);
            w2_g = 0.5 * (x_g + 2 * x * x_g - 3 * x * x * x_g);
            w3_g = 0.5 * x * x * x_g;
            *grad = w0_g * ctrl_pts.row(k - 3) + w1_g * ctrl_pts.row(k - 2) + w2_g * ctrl_pts.row(k - 1) + w3_g * ctrl_pts.row(k);
        }
        if (grad2) {
            *grad2 = Vector4d(w0, w1, w2, w3);
            *idx = k - 3;
        }
        break;
    default:
        throw "getBsplineValueFast parameter error";
        exit(0);
        // ret = Vector3d(INFINITY, INFINITY, INFINITY);
    }

    return ret;
}

inline double UniformBspline::getArcLength(double ts, const MatrixXd &deriv_ctrl_pts, double dt, double start_t, double end_t, vector<double> *arc_lengths) {
    double l = 0.0;
    for (double t = start_t; t < end_t; t += dt) {
        l += getBsplineValueFast<Vector3d>(ts, deriv_ctrl_pts, t, 2).norm() * dt;
        arc_lengths->push_back(l);
    }
    return l;
}