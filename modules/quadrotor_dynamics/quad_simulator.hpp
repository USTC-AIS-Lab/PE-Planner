#pragma once

#include "nominal_quad_dynamic.hpp"

class QuadSimulator {
private:
    Matrix<double, 3, 1> pos_, vel_, acc_, rate_;
    Matrix<double, 4, 1> quat_;
    NominalQuadDynamic quaddynamic_;

public:
    QuadSimulator(const double hover_ratio) : quaddynamic_(hover_ratio) {}
    void set_rate_with_trust(double rx, double ry, double rz , double thrust) {
        Matrix<double, quaddynamic_.x_dim_, 1> state0, state1;
        state0 << pos_, vel_, quat_.w(), quat_.x(), quat_.y(), quat_.z();
        VectorXd u(4);
        u << rx, ry, rz, thrust;
        Matrix<double, NominalQuadDynamic::x_dim_, 1> xdot;
        quaddynamic_.xdot_func(state0, u, Vector3d(0, 0, 0), xdot);
        quaddynamic_.rk4_func(state0, u, Vector3d(0, 0, 0), 0.02, state1);
        pos_ = Matrix<double, 3, 1>(state1[0], state1[1], state1[2]);
        vel_ = Matrix<double, 3, 1>(state1[3], state1[4], state1[5]);
        quat_ = Matrix<double, 4, 1>(state1[7], state1[8], state1[9], state1[6]);
        rate_ = Matrix<double, 3, 1>(rx, ry, rz);
        acc_ = Matrix<double, 3, 1>(xdot[3], xdot[4], xdot[5]);
    }
    void set_pos(double x, double y, double z, double yaw) {
        pos_ = Matrix<double, 3, 1>(x, y, z);
        vel_ = Matrix<double, 3, 1>(0, 0, 0);
        Quaterniond q = AngleAxisd(yaw, Vector3d::UnitZ()) * AngleAxisd(0, Vector3d::UnitY()) * AngleAxisd(0, Vector3d::UnitX());
        quat_ = Matrix<double, 4, 1>(q.x(), q.y(), q.z(), q.w());
    }
    void arm() {

    }
    void disarm() {

    }
    void set_px4_mode(string mode) {

    }
    const Matrix<double, 3, 1> &pos() {return pos_;}
    const Matrix<double, 3, 1> &vel() {return vel_;}
    const Matrix<double, 3, 1> &acc() {return acc_;}
    const Matrix<double, 4, 1> &quat() {return quat_;}
    const Matrix<double, 3, 1> &rate() {return rate_;}
};
