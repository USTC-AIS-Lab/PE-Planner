#include <iostream>

#include "px4_interface.hpp"

using namespace std;
using namespace Eigen;

Px4Interface::Px4Interface(string topic_prefix) {
    ros::NodeHandle nh;
    px4_armed_ = false;
    pose_sub_ = nh.subscribe<geometry_msgs::PoseStamped>(topic_prefix + "/mavros/local_position/pose", 2, boost::bind(&Px4Interface::local_pose_callback, this, _1));
    vel_sub_ = nh.subscribe<geometry_msgs::TwistStamped>(topic_prefix + "/mavros/local_position/velocity_local", 2, boost::bind(&Px4Interface::local_vel_callback, this, _1));
    acc_sub_ = nh.subscribe<sensor_msgs::Imu>(topic_prefix + "/mavros/imu/data_raw", 2, boost::bind(&Px4Interface::acc_callback, this, _1));
    rate_sub_ = nh.subscribe<geometry_msgs::TwistStamped>(topic_prefix + "/mavros/local_position/velocity_body" , 2, boost::bind(&Px4Interface::local_rate_callback, this, _1));
    px4_state_sub_ = nh.subscribe<mavros_msgs::State>(topic_prefix + "/mavros/state", 2, boost::bind(&Px4Interface::px4_state_callback, this, _1));
    pos_target_pub_ = nh.advertise<mavros_msgs::PositionTarget>(topic_prefix + "/mavros/setpoint_raw/local", 1);
    att_target_pub_ = nh.advertise<mavros_msgs::AttitudeTarget>(topic_prefix + "/mavros/setpoint_raw/attitude", 1);
    px4_mode_srv_ = nh.serviceClient<mavros_msgs::SetMode>(topic_prefix + "/mavros/set_mode");
    px4_arm_srv_ = nh.serviceClient<mavros_msgs::CommandBool>(topic_prefix + "/mavros/cmd/arming");
}

void Px4Interface::set_rate_with_trust(double rx, double ry, double rz , double thrust) {
    mavros_msgs::AttitudeTarget cmd;
    cmd.header.stamp = ros::Time::now();
    cmd.header.frame_id = "FCU";
    cmd.body_rate.x = rx;
    cmd.body_rate.y = ry;
    cmd.body_rate.z = rz;
    cmd.thrust = thrust;
    cmd.type_mask = mavros_msgs::AttitudeTarget::IGNORE_ATTITUDE;
    att_target_pub_.publish(cmd);
}

void Px4Interface::set_attitude_with_trust(Quaterniond q, double thrust) {
    mavros_msgs::AttitudeTarget cmd;
    cmd.header.stamp = ros::Time::now();
    cmd.header.frame_id = "FCU";
    cmd.orientation.x = q.x();
	cmd.orientation.y = q.y();
	cmd.orientation.z = q.z();
	cmd.orientation.w = q.w();
    cmd.thrust = thrust;
    cmd.type_mask = mavros_msgs::AttitudeTarget::IGNORE_ROLL_RATE |
					mavros_msgs::AttitudeTarget::IGNORE_PITCH_RATE |
					mavros_msgs::AttitudeTarget::IGNORE_YAW_RATE;
    att_target_pub_.publish(cmd);
}

void Px4Interface::set_pos(double x, double y, double z, double yaw) {
    mavros_msgs::PositionTarget cmd;
    cmd.header.stamp = ros::Time::now();
    cmd.header.frame_id = "FCU";
    cmd.coordinate_frame = 1;
    cmd.position.x = x;
    cmd.position.y = y;
    cmd.position.z = z;
    cmd.yaw = yaw;
    cmd.type_mask = mavros_msgs::PositionTarget::IGNORE_VX + mavros_msgs::PositionTarget::IGNORE_VY + mavros_msgs::PositionTarget::IGNORE_VZ +
                    mavros_msgs::PositionTarget::IGNORE_AFX + mavros_msgs::PositionTarget::IGNORE_AFY + mavros_msgs::PositionTarget::IGNORE_AFZ +
                    mavros_msgs::PositionTarget::FORCE + mavros_msgs::PositionTarget::IGNORE_YAW_RATE;
    pos_target_pub_.publish(cmd);
}

void Px4Interface::arm() {
    mavros_msgs::CommandBool req;
    req.request.value = true;
    if (px4_arm_srv_.call(req)) {
        cout << "Arm" << endl;
    } else {
        cout << "Vehicle arming failed" << endl;
    }
}

void Px4Interface::disarm() {
    mavros_msgs::CommandBool req;
    req.request.value = false;
    if (px4_arm_srv_.call(req)) {
        cout << "Disarm" << endl;
    } else {
        cout << "Vehicle disarming failed" << endl;
    }
}

void Px4Interface::set_px4_mode(string mode) {
    mavros_msgs::SetMode req;
    req.request.custom_mode = mode;
    if (px4_mode_srv_.call(req)) {
        cout << "Set mode: " << mode << endl;
    } else {
        cout << "Failed to set mode: " << mode << endl;
    }
}

void Px4Interface::local_pose_callback(geometry_msgs::PoseStampedConstPtr msg) {
    pos_.x() = msg.get()->pose.position.x;
    pos_.y() = msg.get()->pose.position.y;
    pos_.z() = msg.get()->pose.position.z;
    quat_.x() = msg.get()->pose.orientation.x;
    quat_.y() = msg.get()->pose.orientation.y;
    quat_.z() = msg.get()->pose.orientation.z;
    quat_.w() = msg.get()->pose.orientation.w;
}

void Px4Interface::local_vel_callback(geometry_msgs::TwistStampedConstPtr msg) {
    vel_.x() = msg.get()->twist.linear.x;
    vel_.y() = msg.get()->twist.linear.y;
    vel_.z() = msg.get()->twist.linear.z;
    vel_stamp_ = msg.get()->header.stamp;
}

void Px4Interface::acc_callback(sensor_msgs::ImuConstPtr msg) {
    acc_.x() = msg.get()->linear_acceleration.x;
    acc_.y() = msg.get()->linear_acceleration.y;
    acc_.z() = msg.get()->linear_acceleration.z;
}

void Px4Interface::local_rate_callback(geometry_msgs::TwistStampedConstPtr msg) {
    rate_.x() = msg.get()->twist.angular.x;
    rate_.y() = msg.get()->twist.angular.y;
    rate_.z() = msg.get()->twist.angular.z;
}

void Px4Interface::px4_state_callback(mavros_msgs::StateConstPtr msg) {
    px4_armed_ = msg.get()->armed;
    px4_mode_ = msg.get()->mode;
}
