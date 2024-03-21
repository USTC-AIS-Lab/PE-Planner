#ifndef _PX4_INTERFACE_HPP
#define _PX4_INTERFACE_HPP

#include <string>

#include <Eigen/Dense>
#include <ros/ros.h>

#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <mavros_msgs/PositionTarget.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/CommandBool.h>

using namespace std;
using namespace Eigen;

class Px4Interface {
private:
    //订阅状态信息来源
    ros::Subscriber pose_sub_;
    ros::Subscriber vel_sub_;
    ros::Subscriber acc_sub_;
    ros::Subscriber rate_sub_;
    ros::Subscriber px4_state_sub_;
    //指令话题
    ros::Publisher pos_target_pub_;
    ros::Publisher att_target_pub_;
    ros::ServiceClient px4_mode_srv_;
    ros::ServiceClient px4_arm_srv_;

    Matrix<double, 3, 1> pos_, vel_, acc_, rate_;
    Matrix<double, 4, 1> quat_;
    bool px4_armed_;
    string px4_mode_;
    ros::Time vel_stamp_;

    void local_pose_callback(geometry_msgs::PoseStampedConstPtr msg);
    void local_vel_callback(geometry_msgs::TwistStampedConstPtr msg);
    void acc_callback(sensor_msgs::ImuConstPtr msg);
    void local_rate_callback(geometry_msgs::TwistStampedConstPtr msg);
    void px4_state_callback(mavros_msgs::StateConstPtr msg);

public:
    Px4Interface(string topic_prefix);
    void set_rate_with_trust(double rx, double ry, double rz, double thrust);
    void set_attitude_with_trust(Quaterniond q, double thrust);
    void set_pos(double x, double y, double z, double yaw);
    void arm();
    void disarm();
    void set_px4_mode(string mode);
    const Matrix<double, 3, 1> &pos() {return pos_;}
    const Matrix<double, 3, 1> &vel() {return vel_;}
    const Matrix<double, 3, 1> &acc() {return acc_;}
    const Matrix<double, 4, 1> &quat() {return quat_;}
    const Matrix<double, 3, 1> &rate() {return rate_;}
    ros::Time &vel_stamp() {return vel_stamp_;}
};

#endif