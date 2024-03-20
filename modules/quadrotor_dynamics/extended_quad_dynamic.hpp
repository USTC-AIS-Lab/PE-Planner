#ifndef _NOMINAL_QUAD_DYNAMIC_HPP
#define _NOMINAL_QUAD_DYNAMIC_HPP

#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define PX 0
#define PY 1
#define PZ 2
#define VX 3
#define VY 4
#define VZ 5
#define QW 6
#define QX 7
#define QY 8
#define QZ 9
#define NTH 10
#define RX 0
#define RY 1
#define RZ 2
#define DNTH 3

/*
    A auxiliary dynamic for collective thrust T: T_{k+1}=T_{k}+\Delta T_{k} * \Delta t
    is added to the quadrotor dynamics. 
*/
class ExtendedQuadDynamic {
public:
    static constexpr int x_dim_ = 11; //position, velocity, quaternion and collective thrust (virtual state)
    static constexpr int u_dim_ = 4; //angular, rate of collective thrust (virtual control input)

protected:
    const double hover_ratio_;

    Matrix<double, x_dim_, 1> aux;
    Matrix<double, x_dim_, 1> xdot1;
    Matrix<double, x_dim_, 1> xdot2;
    Matrix<double, x_dim_, 1> xdot3;
    Matrix<double, x_dim_, 1> xdot4;
    Matrix<double, x_dim_, x_dim_> xd1gx0;
    Matrix<double, x_dim_, x_dim_> xd2gx1;
    Matrix<double, x_dim_, x_dim_> xd3gx2;
    Matrix<double, x_dim_, x_dim_> xd4gx3;
    Matrix<double, x_dim_, u_dim_> xd1gu;
    Matrix<double, x_dim_, u_dim_> xd2gu;
    Matrix<double, x_dim_, u_dim_> xd3gu;
    Matrix<double, x_dim_, u_dim_> xd4gu;
    Matrix<double, x_dim_, x_dim_> x1gx0;
    Matrix<double, x_dim_, x_dim_> x2gx0;
    Matrix<double, x_dim_, x_dim_> x3gx0;
    Matrix<double, x_dim_, u_dim_> x1gu;
    Matrix<double, x_dim_, u_dim_> x2gu;
    Matrix<double, x_dim_, u_dim_> x3gu;
    Matrix<double, x_dim_, x_dim_> iden;

public:
    ExtendedQuadDynamic(double hover_ratio) :
        hover_ratio_(hover_ratio) {
        iden = MatrixXd::Identity(x_dim_, x_dim_);
        aux.setZero();
        xdot1.setZero();
        xdot2.setZero();
        xdot3.setZero();
        xdot4.setZero();
        xd1gx0.setZero();
        xd2gx1.setZero();
        xd3gx2.setZero();
        xd4gx3.setZero();
        xd1gu.setZero();
        xd2gu.setZero();
        xd3gu.setZero();
        xd4gu.setZero();
        x1gx0.setZero();
        x2gx0.setZero();
        x3gx0.setZero();
        x1gu.setZero();
        x2gu.setZero();
        x3gu.setZero();
    }

    //xdot = f(x,u), gx = df / dx, gu = df / du
    void xdot_func(const Matrix<double, x_dim_, 1> &x, const VectorXd &u, const Vector3d &disturbance_acc,
        Matrix<double, x_dim_, 1> &xdot, Matrix<double, x_dim_, x_dim_> &gx, 
        Matrix<double, x_dim_, u_dim_> &gu) {
        const double &px = x(0);
        const double &py = x(1);
        const double &pz = x(2);
        const double &vx = x(3);
        const double &vy = x(4);
        const double &vz = x(5);
        const double &qw = x(6);
        const double &qx = x(7);
        const double &qy = x(8);
        const double &qz = x(9);
        const double &nth = x(10); //normalized thrust
        const double &rx = u[0];
        const double &ry = u[1];
        const double &rz = u[2];
        const double &dnth = u[3]; //rate of normalized thrust
        const double thrust = nth * 9.81 / hover_ratio_;

        xdot << vx, vy, vz
                , 2 * (qx * qz + qw * qy) * thrust + disturbance_acc.x()
                , 2 * (qy * qz - qw * qx) * thrust + disturbance_acc.y()
                , ((1 - 2 * (qx * qx + qy * qy)) * thrust - 9.81 + disturbance_acc.z())
                , 0.5 * (-rx * qx - ry * qy - rz * qz)
                , 0.5 * (rx * qw + rz * qy - ry * qz)
                , 0.5 * (ry * qw - rz * qx + rx * qz)
                , 0.5 * (rz * qw + ry * qx - rx * qy)
                , dnth;
        gx(PX, VX) = 1;
        gx(PY, VY) = 1;
        gx(PZ, VZ) = 1;
        gx(VX, QX) = 2 * qz * thrust, gx(VX, QY) = 2 * qw * thrust, gx(VX, QZ) = 2 * qx * thrust, gx(VX, QW) = 2 * qy * thrust;
        gx(VY, QX) = -2 * qw * thrust, gx(VY, QY) = 2 * qz * thrust, gx(VY, QZ) = 2 * qy * thrust, gx(VY, QW) = -2 * qx * thrust;
        gx(VZ, QX) = -4 * qx * thrust, gx(VZ, QY) = -4 * qy * thrust;
        gx(VX, NTH) = 2 * (qx * qz + qw * qy) * 9.81 / hover_ratio_;
        gx(VY, NTH) = 2 * (qy * qz - qw * qx) * 9.81 / hover_ratio_;
        gx(VZ, NTH) = (1 - 2 * (qx * qx + qy * qy)) * 9.81 / hover_ratio_;
        gx(QW, QX) = 0.5 * -rx, gx(QW, QY) = 0.5 * -ry, gx(QW, QZ) = 0.5 * -rz;
        gx(QX, QW) = 0.5 * rx, gx(QX, QY) = 0.5 * rz, gx(QX, QZ) = 0.5 * -ry;
        gx(QY, QW) = 0.5 * ry, gx(QY, QX) = 0.5 * -rz, gx(QY, QZ) = 0.5 * rx;
        gx(QZ, QW) = 0.5 * rz, gx(QZ, QX) = 0.5 * ry, gx(QZ, QY) = 0.5 * -rx;
        // gu(VX, TH) = 2 * (qx * qz + qw * qy) * /*1.084e-5 * 2 * (u[3] * 1000 + 100) * 1000 * 4 / 0.74;//*/9.81 / hover_ratio_;
        // gu(VY, TH) = 2 * (qy * qz - qw * qx) * /*1.084e-5 * 2 * (u[3] * 1000 + 100) * 1000 * 4 / 0.74;//*/9.81 / hover_ratio_;
        // gu(VZ, TH) = (1 - 2 * (qx * qx + qy * qy)) * /*1.084e-5 * 2 * (u[3] * 1000 + 100) * 1000 * 4 / 0.74;//*/9.81 / hover_ratio_;
        gu(NTH, DNTH) = 1;
        gu(QW, RX) = 0.5 * -qx, gu(QW, RY) = 0.5 * -qy, gu(QW, RZ) = 0.5 * -qz;
        gu(QX, RX) = 0.5 * qw, gu(QX, RY) = 0.5 * -qz, gu(QX, RZ) = 0.5 * qy;
        gu(QY, RX) = 0.5 * qz, gu(QY, RY) = 0.5 * qw, gu(QY, RZ) = 0.5 * -qx;
        gu(QZ, RX) = 0.5 * -qy, gu(QZ, RY) = 0.5 * qx, gu(QZ, RZ) = 0.5 * qw;
    }

    void xdot_func(const Matrix<double, x_dim_, 1> &x, const VectorXd &u, const Vector3d &disturbance_acc,
        Matrix<double, x_dim_, 1> &xdot) {
        const double &px = x(0);
        const double &py = x(1);
        const double &pz = x(2);
        const double &vx = x(3);
        const double &vy = x(4);
        const double &vz = x(5);
        const double &qw = x(6);
        const double &qx = x(7);
        const double &qy = x(8);
        const double &qz = x(9);
        const double &nth = x(10); //normalized thrust
        const double &rx = u[0];
        const double &ry = u[1];
        const double &rz = u[2];
        const double &dnth = u[3]; //rate of normalized thrust
        const double thrust = nth * 9.81 / hover_ratio_;
        xdot << vx, vy, vz
                , 2 * (qx * qz + qw * qy) * thrust + disturbance_acc.x()
                , 2 * (qy * qz - qw * qx) * thrust + disturbance_acc.y()
                , ((1 - 2 * (qx * qx + qy * qy)) * thrust - 9.81 + disturbance_acc.z())
                , 0.5 * (-rx * qx - ry * qy - rz * qz)
                , 0.5 * (rx * qw + rz * qy - ry * qz)
                , 0.5 * (ry * qw - rz * qx + rx * qz)
                , 0.5 * (rz * qw + ry * qx - rx * qy)
                , dnth;    
    }

    //x1 = f_rk4(x0, u, dt), gx0 = dx1 / dx0, gu = dx1 / du
    void rk4_func(const Matrix<double, x_dim_, 1> &x0, const VectorXd &u, const Vector3d &disturbance_acc,
        const double &dt, Matrix<double, x_dim_, 1> &x1, 
        Matrix<double, x_dim_, x_dim_> &gx0, Matrix<double, x_dim_, u_dim_> &gu,
        Matrix<double, 3, 1> &acc, Matrix<double, 3, x_dim_> &accdotx0, Matrix<double, 3, u_dim_> &accdotu) {
        xdot_func(x0, u, disturbance_acc, xdot1, xd1gx0, xd1gu);
        acc = xdot1.block(3, 0, 3, 1);
        accdotx0 = xd1gx0.block(3, 0, 3, x_dim_);
        accdotu = xd1gu.block(3, 0, 3, u_dim_);
        aux = x0 + xdot1 * dt / 2.0;
        x1gx0 = iden + dt * 1 / 2.0 * xd1gx0;
        x1gu = dt / 2.0 * xd1gu;
        xdot_func(aux, u, disturbance_acc, xdot2, xd2gx1, xd2gu);
        aux = x0 + xdot2 * dt / 2.0;
        x2gx0 = iden + dt * 1 / 2.0 * xd2gx1 * x1gx0;
        x2gu = dt / 2.0 * (xd2gx1 * x1gu + xd2gu);
        xdot_func(aux, u, disturbance_acc, xdot3, xd3gx2, xd3gu);
        aux = x0 + xdot3 * dt;
        x3gx0 = iden + dt * xd3gx2 * x2gx0;
        x3gu = dt * (xd3gx2 * x2gu + xd3gu);
        xdot_func(aux, u, disturbance_acc, xdot4, xd4gx3, xd4gu);
        x1 = x0 + dt * (1 / 6. * xdot1 + 1 / 3. * xdot2 + 1 / 3. * xdot3 + 1 / 6. * xdot4);
        
        gx0 = iden + 
            dt * (1 / 6.0 * xd1gx0
            + 1 / 3.0 * xd2gx1 * x1gx0
            + 1 / 3.0 * xd3gx2 * x2gx0
            + 1 / 6.0 * xd4gx3 * x3gx0
            );
        gu = dt * (1 / 6.0 * xd1gu + 1 / 3.0 * (xd2gu + xd2gx1 * x1gu) + 1 / 3.0 * (xd3gu + xd3gx2 * x2gu) + 1 / 6.0 * (xd4gu + xd4gx3 * x3gu));
    }

    void rk4_func(const Matrix<double, x_dim_, 1> &x0, const VectorXd &u, const Vector3d &disturbance_acc,
        const double &dt, Matrix<double, x_dim_, 1> &x1) {
        xdot_func(x0, u, disturbance_acc, xdot1);
        aux = x0 + xdot1 * (dt / 2.0);
        xdot_func(aux, u, disturbance_acc, xdot2);
        aux = x0 + xdot2 * (dt / 2.0);
        xdot_func(aux, u, disturbance_acc, xdot3);
        aux = x0 + xdot3 * dt;
        xdot_func(aux, u, disturbance_acc, xdot4);
        x1 = x0 + (dt / 6.0) * (xdot1 + 2 * xdot2 + 2 * xdot3 + xdot4);
    }
};

#endif