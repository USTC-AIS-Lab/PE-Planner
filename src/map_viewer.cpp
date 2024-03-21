#include "map/map.hpp"
#include "ros_interface/ros_interface.hpp"

int main(int argc, char **argv) {
    ros::init(argc, argv, "map_viewer");

    //构造地图
    GridMap gridmap(argv[1]);
    gridmap.update_grid_map();

    //ros接口
    RosInterface ros_inte(gridmap.size());
    sleep(1);
    ros_inte.publish_grid_map(gridmap);
    
    return 0;
}