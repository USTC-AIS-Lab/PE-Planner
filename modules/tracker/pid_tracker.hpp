#pragma once

#include <iostream>
#include <queue>
#include <Eigen/Dense>
#include <ros/ros.h>

using namespace std;
using namespace Eigen;

class PidTracker {
public:
    //param
    const Vector3d p_gain_;//Kp0, Kp1, Kp2
    const Vector3d v_gain_;//Kv0, Kv1, Kv2
    const double rho2_ = 0.998;

    double thr2acc_;
    double P_;
    std::queue<std::pair<ros::Time, double>> timed_thrust_;

    PidTracker(Vector3d p_gain, Vector3d v_gain) :
        p_gain_(p_gain), v_gain_(v_gain) {
        thr2acc_ = 9.81 / 0.309;
        P_ = 1e6;
    }

    double fromQuaternion2yaw(Eigen::Quaterniond q){
        double yaw = atan2(2 * (q.x()*q.y() + q.w()*q.z()), q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z());
        return yaw;
    }

    double computeDesiredCollectiveThrustSignal(const Eigen::Vector3d &des_acc) {
        double throttle_percentage(0.0);
        throttle_percentage = des_acc(2) / thr2acc_;
        return throttle_percentage;
    }

    pair<Quaterniond, double> calculate_control(
        const Vector3d pd, const Vector3d vd, const Vector3d ad,
        Vector3d p, Vector3d v, Quaterniond q) {
        Eigen::Vector3d des_acc(0.0, 0.0, 0.0);
        des_acc = ad + v_gain_.asDiagonal() * (vd - v) + p_gain_.asDiagonal() * (pd - p);
        des_acc += Eigen::Vector3d(0, 0, 9.81);
        double thr = computeDesiredCollectiveThrustSignal(des_acc);
        double roll,pitch;
        double yaw_odom = fromQuaternion2yaw(q);
        double sin = std::sin(yaw_odom);
        double cos = std::cos(yaw_odom);
        roll = (des_acc(0) * sin - des_acc(1) * cos) / 9.81; //小角度假设（arctan(x)=x）
        pitch = (des_acc(0) * cos + des_acc(1) * sin) / 9.81;
        Eigen::Quaterniond qd = Eigen::AngleAxisd(yaw_odom,Eigen::Vector3d::UnitZ())
            * Eigen::AngleAxisd(pitch,Eigen::Vector3d::UnitY())
            * Eigen::AngleAxisd(roll,Eigen::Vector3d::UnitX());
        
        timed_thrust_.push(std::pair<ros::Time, double>(ros::Time::now(), thr));
        while (timed_thrust_.size() > 100)
        {
            timed_thrust_.pop();
        }
        return make_pair(qd, thr);
    }

    bool estimateThrustModel(const Eigen::Vector3d &est_a){
        // cout << "************" << endl;
        ros::Time t_now = ros::Time::now();
        while (timed_thrust_.size() >= 1) {
            // Choose data before 35~45ms ago
            std::pair<ros::Time, double> t_t = timed_thrust_.front();
            double time_passed = (t_now - t_t.first).toSec();
            if (time_passed > 0.045) // 45ms
            {
            // printf("continue, time_passed=%f\n", time_passed);
            timed_thrust_.pop();
            continue;
            }
            if (time_passed < 0.035) // 35ms
            {
            // printf("skip, time_passed=%f\n", time_passed);
            return false;
            }

            // cout << time_passed << endl;

            /***********************************************************/
            /* Recursive least squares algorithm with vanishing memory */
            /***********************************************************/
            double thr = t_t.second;
            timed_thrust_.pop();
            
            /***********************************/
            /* Model: est_a(2) = thr1acc_ * thr */
            /***********************************/
            // cout << thr2acc_ << endl;
            double gamma = 1 / (rho2_ + thr * P_ * thr);
            double K = gamma * P_ * thr;
            double past_thr2acc = thr2acc_;
            thr2acc_ = thr2acc_ + K * (est_a(2) + 9.81 - thr * thr2acc_);
            P_ = (1 - K * thr) * P_ / rho2_;
            // cout << thr2acc_ << " "
            //     << K << " "
            //     << P_ << " "
            //     << gamma << " "
            //     << (est_a(2) + 9.81 - thr * past_thr2acc) << endl;

            // debug_msg_.thr2acc = thr2acc_;
            return true;
        }
        return false;
    }

};
