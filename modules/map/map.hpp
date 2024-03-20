#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class ObsCylinder {
public:
    Vector3d pos_;
    double high_;
    double radius_;

    ObsCylinder(Vector3d pos, double high, double radius) 
        : pos_(pos), high_(high), radius_(radius) {}
};

class ObsRing {
public:
    Vector3d pos_;
    double wid_; //width of ring
    double inner_r_; //inner radius

    ObsRing(Vector3d pos, double wid, double inner_r) 
        : pos_(pos), wid_(wid), inner_r_(inner_r) {}
};

class ObsWall {
public:
    Vector3d p1_;
    Vector3d p2_;

    ObsWall(Vector3d p1, Vector3d p2)
        : p1_(p1), p2_(p2) {}
};

class GridMap {
private:
    vector<bool> grid_map_;
    double resolution_; //resolution of gridmap
    Vector3d size_; //map size
    Vector3i isize_; //gridmap size

    vector<ObsCylinder> cylinders_;
    vector<ObsRing> rings_;
    vector<ObsWall> walls_;

public:
    GridMap() {
        resolution_ = 0.0;
        size_ = Vector3d(0., 0., 0.);
        isize_ = Vector3i(0., 0., 0.);
    }
    GridMap(double resolution, Vector3d size) : resolution_(resolution), size_(size) {
        isize_ = Vector3i(round(size_(0) / resolution_), round(size_(1) / resolution_), round(size_(2) / resolution_));
        grid_map_.resize(isize_(0) * isize_(1) * isize_(2));
        for (int i = 0; i < grid_map_.size(); i++) {
            grid_map_[i] = false;
        }
    }
    GridMap(string file) {
        read_map_file(file);
    }
    void update_grid_map();
    bool read_map_file(string path);
    Vector3d idx2pos(const Vector3i idx) const {return Vector3d((idx.x() + 0.5) * resolution_, (idx.y() + 0.5) * resolution_, (idx.z() + 0.5) * resolution_);}
    Vector3i pos2idx(const Vector3d pos) const {return Vector3i(pos.x() / resolution_, pos.y() / resolution_, pos.z() / resolution_);}
    bool check_in_map(Vector3i idx) const {
        if (idx.x() >= isize_(0) || idx.y() >= isize_(1) || idx.z() >= isize_(2) 
            || idx.x() < 0 || idx.y() < 0 || idx.z() < 0) {
            return false;
        }
        return true;
    }
    bool check_in_map(Vector3d pos) const {
        Vector3i idx = pos2idx(pos);
        if (idx.x() >= isize_(0) || idx.y() >= isize_(1) || idx.z() >= isize_(2) 
            || idx.x() < 0 || idx.y() < 0 || idx.z() < 0) {
            return false;
        }
        return true;
    }
    bool operator()(Vector3i idx) const {
        if (check_in_map(idx)) {
            return grid_map_[idx.z() * isize_(1) * isize_(0) + idx.y() * isize_(0) + idx.x()];
        } else {
            return true;
        }
    }
    bool operator()(int x, int y, int z) const {
        Vector3i idx(x, y, z);
        if (check_in_map(idx)) {
            return grid_map_[idx.z() * isize_(1) * isize_(0) + idx.y() * isize_(0) + idx.x()];
        } else {
            return true;
        }
    }
    bool operator()(Vector3d pos) const {
        Vector3i idx = pos2idx(pos);
        if (check_in_map(idx)) {
            return grid_map_[idx.z() * isize_(1) * isize_(0) + idx.y() * isize_(0) + idx.x()];
        } else {
            return true;
        }
    }
    bool operator()(double x, double y, double z) const {
        Vector3i idx = pos2idx(Vector3d(x, y, z));
        if (check_in_map(idx)) {
            return grid_map_[idx.z() * isize_(1) * isize_(0) + idx.y() * isize_(0) + idx.x()];
        } else {
            return true;
        }
    }
    void set_grid(Vector3i idx, bool obs) {
        if (check_in_map(idx)) {
            grid_map_[idx.z() * isize_(1) * isize_(0) + idx.y() * isize_(0) + idx.x()] = obs;
        }
    } 
    void set_grid(Vector3d pos, bool obs) {
        Vector3i idx = pos2idx(pos);
        set_grid(idx, obs);
    }
    const double resolution() const {return resolution_;}
    const Vector3d size() const {return size_;}
    const Vector3i isize() const {return isize_;}
    void push_obs(ObsCylinder &o) {
        cylinders_.push_back(o);
    }
    void push_obs(ObsRing &o) {
        rings_.push_back(o);
    }
    void push_obs(ObsWall &o) {
        walls_.push_back(o);
    }
};
