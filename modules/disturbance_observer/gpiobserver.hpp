#pragma once

#include <Eigen/Core>

using namespace Eigen;

class GPIObserver {
private:
    Vector3d v_hat_;
    Vector3d d_hat_;
    Vector3d dd_hat_;
    DiagonalMatrix<double, 3> l1_, l2_, l3_;

public:
    GPIObserver(Vector3d l1, Vector3d l2, Vector3d l3) : l1_(l1), l2_(l2), l3_(l3) {
        v_hat_.setZero();
        d_hat_.setZero();
        dd_hat_.setZero();
    }
    void update(const Vector3d v_measure, const Matrix<double, 3, 3> &r, const double thrust, const double dt) {
        d_hat_ += dd_hat_ * dt;
        v_hat_ += (r * Vector3d(0, 0, thrust) - Vector3d(0, 0, 9.81) + d_hat_) * dt;
        dd_hat_ += l3_ * (v_measure - v_hat_) * dt;
        d_hat_ += l2_ * (v_measure - v_hat_) * dt;
        v_hat_ += l1_ * (v_measure - v_hat_) * dt;
        // Vector3d err = v_measure - v_hat_;
        // v_hat_ += dt * (r * Vector3d(0, 0, thrust) - Vector3d(0, 0, 9.81) + d_hat_ + l1_ * err);
        // d_hat_ += dt * (dd_hat_ + l2_ * err);
        // dd_hat_ += dt * (l3_ * err);
    }
    Vector3d get_dhat() const {
        return d_hat_;
        // return Vector3d(0., 0., 0.);
    }
    Vector3d get_vhat() const {
        return v_hat_;
    }
};