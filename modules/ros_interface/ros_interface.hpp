#pragma once

#include <ros/ros.h>
#include <tf/transform_broadcaster.h>

#include "map/map.hpp"
#include "map/sdf.hpp"

class RosInterface {
private:
    ros::Publisher grid_map_pub_;
    ros::Publisher sdf_map_pub_;
    ros::Publisher quadmesh_pub_;
    ros::Publisher kino_traj_pub_;
    ros::Publisher bspline_traj_pub_;
    ros::Publisher mpcc_traj_pub_;
    ros::Publisher predict_traj_pub_;
    ros::Publisher collision_pub_;
    ros::Publisher fanmesh_pub_;
    ros::Publisher dyn_obs_pub_;
    tf::TransformBroadcaster broadcaster;
    Vector3d map_size_;
public:
    RosInterface(Vector3d map_size);
    void publish_grid_map(GridMap &map);
    void publish_sdf_map(SdfMap &sdf);
    void publish_quadmesh(Vector3d pos, Vector4d quat);
    void publish_kino_traj(vector<Vector3d> &traj);
    void publish_bspline_traj(vector<Vector3d> &traj);
    void publish_mpcc_traj(vector<Vector3d> &traj);
    void publish_predict_traj(vector<Vector3d> &traj);
    void publish_collision(vector<Vector3d> &pos);
    void publish_fanmesh(Vector3d pos, Vector3d ang);
    void publish_dyn_obs(vector<DynObs> &obs);
    void publish_pose(Vector3d pos, Vector4d quat);
};
