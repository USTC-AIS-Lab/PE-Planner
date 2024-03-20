#include <fstream>
#include <stdio.h>

#include "map.hpp"

bool GridMap::read_map_file(string path) {
    ifstream file;
    file.open(path);
    if (!file) {
        cerr << "Fail to open map file" << endl;
        throw "Fail to open map file";
        return false;
    }
    cylinders_.clear();
    rings_.clear();
    walls_.clear();
    string s;
    while (getline(file, s)) {
        char keyc[25];
        sscanf(s.c_str(), "%s", keyc);
        string key(keyc);
        if (key == "map") {
            double v1, v2, v3, v4;
            sscanf(s.c_str() + key.size(), "%lf %lf %lf %lf", &v1, &v2, &v3, &v4);
            resolution_ = v1;
            size_ = Vector3d(v2, v3, v4);
            isize_ = Vector3i(round(size_(0) / resolution_), round(size_(1) / resolution_), round(size_(2) / resolution_));
            grid_map_.resize(isize_(0) * isize_(1) * isize_(2));
            for (int i = 0; i < grid_map_.size(); i++) {
                grid_map_[i] = false;
            }
        } else if (key == "ObsWall") {
            double v1, v2, v3, v4, v5, v6;
            sscanf(s.c_str() + key.size(), "%lf %lf %lf %lf %lf %lf", &v1, &v2, &v3, &v4, &v5, &v6);
            walls_.push_back(ObsWall(Vector3d(v1, v2, v3), Vector3d(v4, v5, v6)));
        } else if (key == "ObsCylinder") {
            double v1, v2, v3, v4, v5;
            sscanf(s.c_str() + key.size(), "%lf %lf %lf %lf %lf", &v1, &v2, &v3, &v4, &v5);
            cylinders_.push_back(ObsCylinder(Vector3d(v1, v2, v3), v4, v5));
        } else if (key == "ObsRing") {
            double v1, v2, v3, v4, v5;
            sscanf(s.c_str() + key.size(), "%lf %lf %lf %lf %lf", &v1, &v2, &v3, &v4, &v5);
            rings_.push_back(ObsRing(Vector3d(v1, v2, v3), v4, v5));
        } else {
            throw "Fail to recognize " + key;
            return false;
        }
    }

    return true;
}

void GridMap::update_grid_map() {
    grid_map_.resize(isize_(0) * isize_(1) * isize_(2));
    for (int i = 0; i < grid_map_.size(); i++) {
        grid_map_[i] = false;
    }
    for (int i = 0; i < cylinders_.size(); i++) {
        const auto &cyl = cylinders_[i];
        for (double r = resolution_ / 2.0; r < cyl.radius_; r += resolution_ / 2.0) {
            for (double ang = 0; ang < M_PI * 2; ang += M_PI / 100) {
                for (double h = 0.0; h < cyl.high_; h += resolution_ / 2.0) {
                    Vector3d pos = cyl.pos_;
                    pos(0) += r * sin(ang);
                    pos(1) += r * cos(ang);
                    pos(2) += h;
                    set_grid(pos, true);
                }
            }
        }
    }
    for (int i = 0; i < walls_.size(); i++) {
        const auto &wall = walls_[i];
        for (double x = wall.p1_(0) + 1e-6; x < wall.p2_(0) - 1e-6; x += resolution_ / 2.0) {
            for (double y = wall.p1_(1) + 1e-6; y < wall.p2_(1) - 1e-6; y += resolution_ / 2.0) {
                for (double z = wall.p1_(2) + 1e-6; z < wall.p2_(2) - 1e-6; z += resolution_ / 2.0) {
                    set_grid(Vector3d(x, y, z), true);
                }
            }
        }
    }
    for (int i = 0; i < rings_.size(); i++) {
        const auto &ring = rings_[i];
        for (double r = 0.0; r < ring.inner_r_; r += resolution_ / 2.0) {
            for (double ang = 0; ang < M_PI * 2; ang += M_PI / 200) {
                for (double w = 0; w < ring.wid_; w += resolution_ / 2.0) {
                    Vector3d pos = ring.pos_;
                    pos(0) += r * sin(ang);
                    pos(2) += r * cos(ang);
                    pos(1) += w - ring.wid_ / 2.0;
                    set_grid(pos, false);
                }
            }
        }
    }
    for (int z = 0; z < isize_(2); z++) {
        for (int x = 0; x < isize_(0); x++) {
            for (int y = 0; y < isize_(1); y++) {
                if (x == 0 || x == isize_(0) - 1 || y == 0 || y == isize_(1) - 1 || z == 0 || z == isize_(2) - 1) {
                    set_grid(Vector3i(x, y, z), true);
                }
            }
        }
    }
}