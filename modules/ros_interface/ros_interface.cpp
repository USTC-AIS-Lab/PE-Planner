#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/Marker.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <tf2/LinearMath/Quaternion.h>

#include "ros_interface.hpp"
#include "common/rotation_math.hpp"

#define COLOR(color, R, G, B, A) color.a = A / 255.0, color.r = R / 255.0, color.g = G / 255.0, color.b = B / 255.0

RosInterface::RosInterface(Vector3d map_size)
    : map_size_(map_size) {
    ros::NodeHandle nh;
    grid_map_pub_ = nh.advertise<sensor_msgs::PointCloud2>("/map", 1);
    sdf_map_pub_ = nh.advertise<sensor_msgs::PointCloud2>("/sdf_map", 1);
    quadmesh_pub_ = nh.advertise<visualization_msgs::Marker>("/quadrotor", 1);
    kino_traj_pub_ = nh.advertise<visualization_msgs::Marker>("/kino_traj", 1);
    bspline_traj_pub_ = nh.advertise<visualization_msgs::Marker>("/bspline_traj", 1);
    mpcc_traj_pub_ = nh.advertise<visualization_msgs::Marker>("/mpcc_traj", 1);
    predict_traj_pub_ = nh.advertise<visualization_msgs::Marker>("/predict_traj", 1);
    collision_pub_ = nh.advertise<visualization_msgs::Marker>("/collision", 1);
    fanmesh_pub_ = nh.advertise<visualization_msgs::Marker>("/fan", 1);
    dyn_obs_pub_ = nh.advertise<sensor_msgs::PointCloud2>("/dyn_obs", 1);
}

void RosInterface::publish_grid_map(GridMap &map) {
    pcl::PointCloud<pcl::PointXYZ> cloudMap;
    double resolution = map.resolution();
    for (double z = resolution / 2.0; z < map_size_(2); z += resolution) {
        for (double y = resolution / 2.0; y < map_size_(1); y += resolution) {
            for (double x = resolution / 2.0; x < map_size_(0); x += resolution) {
                if (map(Vector3d(x, y ,z))
                    && map.pos2idx(Vector3d(x, y, z)).z() != map.isize().z() - 1
                    && map.pos2idx(Vector3d(x, y, z)).z() != 0) {
                    pcl::PointXYZ pt;
                    pt.x = x - map_size_(0) / 2.0, pt.y = y - map_size_(1) / 2.0, pt.z = z;
                    cloudMap.points.push_back(pt);
                }
            }
        }
    }
    cloudMap.width = cloudMap.points.size();
    cloudMap.height = 1;
    cloudMap.is_dense = true;

    sensor_msgs::PointCloud2 globalMap_pcd;
    pcl::toROSMsg(cloudMap, globalMap_pcd);
    globalMap_pcd.header.frame_id = "world";
    grid_map_pub_.publish(globalMap_pcd);
}

void RosInterface::publish_sdf_map(SdfMap &sdf) {
    pcl::PointCloud<pcl::PointXYZI> cloudMap;
    double resolution = sdf.resolution_;
    double z = 1.0;
    // for (double z = resolution / 2.0; z < sdf.map_size_(2); z += resolution) {
        // double y = 1.0;
        for (double y = resolution / 2.0; y < sdf.map_size_(1); y += resolution) {
            for (double x = resolution / 2.0; x < sdf.map_size_(0); x += resolution) {
                pcl::PointXYZI pt;
                pt.x = x - map_size_(0) / 2.0, pt.y = y - map_size_(1) / 2.0, pt.z = 0.0;
                pt.intensity = sdf.get_dist_with_grad_trilinear(Vector3d(x, y, z)).first;
                cloudMap.points.push_back(pt);
            }
        }
    // }
    cloudMap.width = cloudMap.points.size();
    cloudMap.height = 1;
    cloudMap.is_dense = true;

    sensor_msgs::PointCloud2 globalMap_pcd;
    pcl::toROSMsg(cloudMap, globalMap_pcd);
    globalMap_pcd.header.frame_id = "world";
    sdf_map_pub_.publish(globalMap_pcd);
}

void RosInterface::publish_quadmesh(Vector3d pos, Vector4d quat) {
    visualization_msgs::Marker meshROS;
    meshROS.header.frame_id = "world";
    meshROS.header.stamp = ros::Time::now(); 
    meshROS.ns = "quadmesh";
    meshROS.id = 0;
    meshROS.type = visualization_msgs::Marker::MESH_RESOURCE;
    meshROS.action = visualization_msgs::Marker::ADD;
    meshROS.pose.position.x = pos.x() - map_size_(0) / 2.0;
    meshROS.pose.position.y = pos.y() - map_size_(1) / 2.0;
    meshROS.pose.position.z = pos.z();
    Vector3d ang = Quaterniond(quat.w(), quat.x(), quat.y(), quat.z()).matrix().eulerAngles(2, 1, 0);
    ang(0) += M_PI / 4;
    Quaterniond q = AngleAxisd(ang[0], Vector3d::UnitZ()) * AngleAxisd(ang[1], Vector3d::UnitY()) * AngleAxisd(ang[2], Vector3d::UnitX());
    meshROS.pose.orientation.w = q.w();
    meshROS.pose.orientation.x = q.x();
    meshROS.pose.orientation.y = q.y();
    meshROS.pose.orientation.z = q.z();
    meshROS.scale.x = 1.0;
    meshROS.scale.y = 1.0;
    meshROS.scale.z = 1.0;
    // meshROS.color.a = 1.0;
    // meshROS.color.r = 1.0;
    // meshROS.color.g = 0.0;
    // meshROS.color.b = 0.0;
    COLOR(meshROS.color, 255, 20, 147, 255);
    meshROS.mesh_resource = std::string("file:///home/jiaxin/workspace/PX4/quadrotor_mpc/meshes/hummingbird.mesh");
    quadmesh_pub_.publish(meshROS);     
}

void RosInterface::publish_fanmesh(Vector3d pos, Vector3d ang) {
    visualization_msgs::Marker meshROS;
    meshROS.header.frame_id = "world";
    meshROS.header.stamp = ros::Time::now(); 
    meshROS.ns = "fenmesh";
    meshROS.id = 0;
    meshROS.type = visualization_msgs::Marker::MESH_RESOURCE;
    meshROS.action = visualization_msgs::Marker::ADD;
    meshROS.pose.position.x = pos.x() - map_size_(0) / 2.0;
    meshROS.pose.position.y = pos.y() - map_size_(1) / 2.0;
    meshROS.pose.position.z = pos.z();
    Quaterniond q = AngleAxisd(ang[2], Vector3d::UnitZ()) * AngleAxisd(ang[1], Vector3d::UnitY()) * AngleAxisd(ang[0], Vector3d::UnitX());
    meshROS.pose.orientation.w = q.w();
    meshROS.pose.orientation.x = q.x();
    meshROS.pose.orientation.y = q.y();
    meshROS.pose.orientation.z = q.z();
    meshROS.scale.x = 0.002;
    meshROS.scale.y = 0.002;
    meshROS.scale.z = 0.002;
    meshROS.color.a = 1.0;
    meshROS.color.r = 0.0;
    meshROS.color.g = 1.0;
    meshROS.color.b = 0.0;
    meshROS.mesh_resource = std::string("file:///home/jiaxin/workspace/PX4/quadrotor_mpc/meshes/fan.stl");
    fanmesh_pub_.publish(meshROS);     
}

void RosInterface::publish_kino_traj(vector<Vector3d> &traj) {
    visualization_msgs::Marker tr;
    tr.header.frame_id = "world";
    tr.header.stamp = ros::Time::now(); 
    tr.ns = "kino trajectory";
    tr.action = visualization_msgs::Marker::ADD;
    tr.pose.orientation.w = 1.0;
    tr.pose.orientation.x = 0.0;
    tr.pose.orientation.y = 0.0;
    tr.pose.orientation.z = 0.0;
    tr.id    = 0;
    tr.type    = visualization_msgs::Marker::LINE_STRIP;
    tr.scale.x = 0.05 * 1.5;
    tr.scale.y = 0.05 * 1.5;
    tr.scale.z = 0.05 * 1.5;
    COLOR(tr.color, 0, 203, 0, 200);
    // tr.color.a = 1.0;
    // tr.color.r = 0.0;
    // tr.color.g = 1.0;
    // tr.color.b = 0.0;
    for (const auto &tp : traj) {
        geometry_msgs::Point p;
        p.x = tp.x() - map_size_(0) / 2.0, p.y = tp.y() - map_size_(1) / 2.0, p.z = tp.z();
        tr.points.push_back(p);
    }
    kino_traj_pub_.publish(tr);
}

void RosInterface::publish_bspline_traj(vector<Vector3d> &traj) {
    visualization_msgs::Marker tr;
    tr.header.frame_id = "world";
    tr.header.stamp = ros::Time::now(); 
    tr.ns = "bspline trajectory";
    tr.action = visualization_msgs::Marker::ADD;
    tr.pose.orientation.w = 1.0;
    tr.pose.orientation.x = 0.0;
    tr.pose.orientation.y = 0.0;
    tr.pose.orientation.z = 0.0;
    tr.id    = 0;
    tr.type    = visualization_msgs::Marker::LINE_STRIP;
    tr.scale.x = 0.05 * 1.5;
    tr.scale.y = 0.05 * 1.5;
    tr.scale.z = 0.05 * 1.5;
    COLOR(tr.color, 205, 133, 63, 255);
    // tr.color.a = 200 / 255.0;
    // tr.color.r = 1.0;
    // tr.color.g = 0.0;
    // tr.color.b = 0.0;
    for (const auto &tp : traj) {
        geometry_msgs::Point p;
        p.x = tp.x() - map_size_(0) / 2.0, p.y = tp.y() - map_size_(1) / 2.0, p.z = tp.z();
        tr.points.push_back(p);
    }
    bspline_traj_pub_.publish(tr);
}

void RosInterface::publish_mpcc_traj(vector<Vector3d> &traj) {
    visualization_msgs::Marker tr;
    tr.header.frame_id = "world";
    tr.header.stamp = ros::Time::now(); 
    tr.ns = "mpcc trajectory";
    tr.action = visualization_msgs::Marker::ADD;
    tr.pose.orientation.w = 1.0;
    tr.pose.orientation.x = 0.0;
    tr.pose.orientation.y = 0.0;
    tr.pose.orientation.z = 0.0;
    tr.id    = 0;
    tr.type    = visualization_msgs::Marker::LINE_STRIP;
    tr.scale.x = 0.05 * 1.5;
    tr.scale.y = 0.05 * 1.5;
    tr.scale.z = 0.05 * 1.5;
    COLOR(tr.color, 30, 144, 255, 255);
    // tr.color.a = 1.0;
    // tr.color.r = 1.0;
    // tr.color.g = 106. / 255.0;
    // tr.color.b = 106. / 255.0;
    for (const auto &tp : traj) {
        geometry_msgs::Point p;
        p.x = tp.x() - map_size_(0) / 2.0, p.y = tp.y() - map_size_(1) / 2.0, p.z = tp.z();
        tr.points.push_back(p);
    }
    mpcc_traj_pub_.publish(tr);
}

void RosInterface::publish_predict_traj(vector<Vector3d> &traj) {
    visualization_msgs::Marker tr;
    tr.header.frame_id = "world";
    tr.header.stamp = ros::Time::now(); 
    tr.ns = "predict trajectory";
    tr.action = visualization_msgs::Marker::ADD;
    tr.pose.orientation.w = 1.0;
    tr.pose.orientation.x = 0.0;
    tr.pose.orientation.y = 0.0;
    tr.pose.orientation.z = 0.0;
    tr.id    = 0;
    tr.type    = visualization_msgs::Marker::LINE_STRIP;
    tr.scale.x = 0.05 * 1.5;
    tr.scale.y = 0.05 * 1.5;
    tr.scale.z = 0.05 * 1.5;
    COLOR(tr.color, 255, 69, 0, 255);
    // tr.color.a = 255 / 255.0;
    // tr.color.r = 138. / 255.;
    // tr.color.g = 43. / 255.;
    // tr.color.b = 226. / 255.;
    for (const auto &tp : traj) {
        geometry_msgs::Point p;
        p.x = tp.x() - map_size_(0) / 2.0, p.y = tp.y() - map_size_(1) / 2.0, p.z = tp.z();
        tr.points.push_back(p);
    }
    predict_traj_pub_.publish(tr);
}

void RosInterface::publish_collision(vector<Vector3d> &pos) {
    visualization_msgs::Marker cp;
    cp.header.frame_id = "world";
    cp.header.stamp = ros::Time::now(); 
    cp.ns = "collision";
    cp.action = visualization_msgs::Marker::ADD;
    cp.pose.orientation.w = 1.0;
    cp.pose.orientation.x = 0.0;
    cp.pose.orientation.y = 0.0;
    cp.pose.orientation.z = 0.0;
    cp.id    = 0;
    cp.type    = visualization_msgs::Marker::SPHERE_LIST;
    cp.scale.x = 0.1 * 2;
    cp.scale.y = 0.1 * 2;
    cp.scale.z = 0.1 * 2;
    cp.color.a = 1.0;
    cp.color.r = 205. / 255.;
    cp.color.g = 133. / 255.;
    cp.color.b = 0. / 255.;
    for (const auto &tp : pos) {
        geometry_msgs::Point p;
        p.x = tp.x() - map_size_(0) / 2.0, p.y = tp.y() - map_size_(1) / 2.0, p.z = tp.z();
        cp.points.push_back(p);
    }
    predict_traj_pub_.publish(cp);
}

void RosInterface::publish_dyn_obs(vector<DynObs> &obs) {
    pcl::PointCloud<pcl::PointXYZ> cloudMap;
    double resolution = 0.05;
    for (const auto &o : obs) {
        for (double z = resolution / 2.0; z < map_size_(2); z += resolution) {
            for (double a = 0; a < 360.0; a += 10.0) {
                for (double r = 0; r < o.radius_; r += resolution / 2.0) {
                    pcl::PointXYZ pt;
                    pt.x = o.pos_.x() + r * sin(a / 180.0 * M_PI) - map_size_(0) / 2.0;
                    pt.y = o.pos_.y() + r * cos(a / 180.0 * M_PI) - map_size_(1) / 2.0;
                    pt.z = z;
                    cloudMap.points.push_back(pt);
                }
            }
        }
    }
    cloudMap.width = cloudMap.points.size();
    cloudMap.height = 1;
    cloudMap.is_dense = true;

    sensor_msgs::PointCloud2 dynobs_pcd;
    pcl::toROSMsg(cloudMap, dynobs_pcd);
    dynobs_pcd.header.frame_id = "world";
    dyn_obs_pub_.publish(dynobs_pcd);
}

void RosInterface::publish_pose(Vector3d pos, Vector4d quat) {
    geometry_msgs::PoseStamped pose;
    pose.header.frame_id = "world";
    pose.header.stamp = ros::Time::now();
    pose.pose.position.x = pos.x() - map_size_(0) / 2.0;
    pose.pose.position.y = pos.y() - map_size_(1) / 2.0;
    pose.pose.position.z = pos.z();
    pose.pose.orientation.x = quat.x();
    pose.pose.orientation.y = quat.y();
    pose.pose.orientation.z = quat.z();
    pose.pose.orientation.w = quat.w();
    tf::Transform transform;
    transform.setOrigin(tf::Vector3(pos.x() - map_size_(0) / 2.0, pos.y() - map_size_(1) / 2.0, pos.z()));
    Vector3d ang = quaternion_to_rpy(Quaterniond(quat.w(), quat.x(), quat.y(), quat.z()));
    ang(1) = ang(0) = 0;
    Quaterniond q = rpy_to_quaternion(ang);
    transform.setRotation(tf::Quaternion(quat.x(), quat.y(), quat.z(), quat.w()));
    broadcaster.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", "quad"));
}
