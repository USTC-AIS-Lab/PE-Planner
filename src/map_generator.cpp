#include <random>
#include <fstream>

#include "map/map.hpp"
#include "ros_interface/ros_interface.hpp"

using namespace std;

int main(int argc, char **argv) {
    ros::init(argc, argv, "map_generator");

    if (argc != 2) {
        return 0;
    }

    const Vector3d map_size(50.0, 10.0, 3.0);
    int obs_num = 50;
    double max_radius = 0.3;
    double min_radius = 0.2;

    //构造地图
    GridMap gridmap(0.05, map_size);

    random_device rd;
    default_random_engine eng(rd());
    auto rand_x = uniform_real_distribution<double>(0, map_size(0));
    auto rand_y = uniform_real_distribution<double>(0, map_size(1));
    auto rand_r = uniform_real_distribution<double>(min_radius, max_radius);
    vector<ObsCylinder> cylinders;
    for(int i = 0; i < obs_num; i++) {
        while (true) {
            double x = rand_x(eng);
            double y = rand_y(eng);
            double r = rand_r(eng);

            // cout << x << " " << y << " " << r << endl;
            if (x + r > map_size(0) - 2 || x - r < 2) {
                continue;
            }
            if (y + r > map_size(1) - 0.9 || y - r < 0.9) {
                continue;
            }
            bool flag = false;
            for (auto &c : cylinders) {
                double dis = (Vector3d(x, y, 0) - c.pos_).norm() - r - c.radius_;
                if (dis < 0.8) {
                    flag = true;
                    break;
                }
            }
            if (flag) {
                continue;
            }

            cylinders.push_back(ObsCylinder(Vector3d(x, y, 0), map_size(2), r));
            break;
        }
    }
    for (auto &c : cylinders) {
        gridmap.push_obs(c);
    }

    gridmap.update_grid_map();

    //构建SDF地图
    SdfMap sdfmap(gridmap.resolution() / 1.0, gridmap);

    //ros接口
    RosInterface ros_inte(gridmap.size());
    sleep(1);
    ros_inte.publish_grid_map(gridmap);
    ros_inte.publish_sdf_map(sdfmap);
    sleep(1);

    ofstream map_file(argv[1]);
    map_file << "map " << 0.05 << " " << map_size(0) << " " << map_size(1) << " " << map_size(2) << endl;
    for (auto &c : cylinders) {
        map_file << "ObsCylinder " << c.pos_(0) << " " << c.pos_(1) << " " << c.pos_(2) << " " << c.high_ << " " << c.radius_ << endl;
    }
    map_file.close();
    
    return 0;
}