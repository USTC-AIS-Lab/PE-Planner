#pragma once

#include <iostream>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

inline Vector3d quaternion_to_rpy(const Quaterniond &q) {
    const double &qw = q.w();
    const double &qx = q.x();
    const double &qy = q.y();
    const double &qz = q.z();

    Vector3d rpy;
    /* roll  */ rpy.x() = std::atan2(2. * (qw*qx + qy*qz), 1. - 2. * (qx*qx + qy*qy));
    double sin_pitch = 2. * (qw*qy - qz*qx);
    // sin_pitch = sin_pitch >  1.0 ?  1.0 : sin_pitch;
    // sin_pitch = sin_pitch < -1.0 ? -1.0 : sin_pitch;
    /* pitch */ rpy.y() = std::asin(sin_pitch);
    /* yaw   */ rpy.z() = std::atan2(2. * (qw*qz + qx*qy), 1. - 2. * (qy*qy + qz*qz));
    return rpy;
}

inline Matrix<double, 3, 3> quaternion_to_matrix(const Quaterniond &q) {
    auto &x = q.x();
    auto &y = q.y();
    auto &z = q.z();
    auto &w = q.w();
    Matrix<double, 3, 3> m;
    m << 1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
        2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
        2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y);
    return m;
}


inline Quaterniond rpy_to_quaternion(Vector3d rpy) // yaw (Z), pitch (Y), roll (X)
{
    // Abbreviations for the various angular functions
    double cy = cos(rpy[2] * 0.5);
    double sy = sin(rpy[2] * 0.5);
    double cp = cos(rpy[1] * 0.5);
    double sp = sin(rpy[1] * 0.5);
    double cr = cos(rpy[0] * 0.5);
    double sr = sin(rpy[0] * 0.5);
 
    Quaterniond q;
    q.w() = cy * cp * cr + sy * sp * sr;
    q.x() = cy * cp * sr - sy * sp * cr;
    q.y() = sy * cp * sr + cy * sp * cr;
    q.z() = sy * cp * cr - cy * sp * sr;
 
    return q;
}