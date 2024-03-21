#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <random>
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

class DynObs {
public:
    double radius_;
    double vel_;
    Vector2d start_pos_;
    Vector2d end_pos_;

    Vector2d pos_;
    double pos_ratio_;
    double vel_ratio_;
    int dir_;

    Vector3d axis_l_;

    DynObs(double radius, double vel, Vector2d start_pos, Vector2d end_pos) 
        : radius_(radius), vel_(vel), start_pos_(start_pos), end_pos_(end_pos) {
        vel_ratio_ = vel_ / (end_pos_ - start_pos_).norm();
        random_device rd;
        default_random_engine eng(rd());
        pos_ratio_ = uniform_real_distribution<double>(0, 1)(eng);
        pos_ = pos_ratio_ * start_pos_ + (1 - pos_ratio_) * end_pos_;
        dir_ = uniform_int_distribution<int>(0, 1)(eng) * 2 - 1;

        double alpha = 4 * radius_ * radius_;
        double beta = 3.0 * 3.0;
        double tmp1 = 2 / alpha;
        double tmp2 = (4.0 - alpha * tmp1) / beta;
        axis_l_ = Vector3d(sqrt(1 / tmp1), sqrt(1 / tmp1), sqrt(1 / tmp2));
    }
    
    void update(double dt) {
        pos_ratio_ += dir_ * vel_ratio_ * dt;
        if (pos_ratio_ < -1e-6) {
            dir_ = 1;
            pos_ratio_ = 0;
        } else if (pos_ratio_ > 1 + 1e-6) {
            dir_ = -1;
            pos_ratio_ = 1;
        }
        pos_ = pos_ratio_ * start_pos_ + (1 - pos_ratio_) * end_pos_;
    }

    void get_dis(const Vector3d pos, const double t, double *dis, Vector3d *g = nullptr) const {
        double r = pos_ratio_ + dir_ * vel_ratio_ * t;
        Vector2d p = r * start_pos_ + (1 - r) * end_pos_;
        double tmp = (pos - Vector3d(p.x(), p.y(), pos.z())).norm();
        *dis = tmp - radius_;
        if (g) {
            if (tmp > 1e-12) {
                *g = (pos - Vector3d(p.x(), p.y(), pos.z())) / tmp;
            } else {
                g->setZero();
            }
        }
    }

    void get_dis_ellipsoid(const Vector3d pos, const double t, double *dis, Vector3d *g = nullptr) const {
        double r = pos_ratio_ + dir_ * vel_ratio_ * t;
        Vector2d _p = r * start_pos_ + (1 - r) * end_pos_;
        Vector3d p(_p.x(), _p.y(), 1.5);
        Vector3d m = pos - p; //ignore rotation
        Vector3d mz(0, 0, m.z());
        Vector3d mxy = m - mz;
        double sin2_psi = mxy.squaredNorm() / m.squaredNorm();
        double cos2_psi = mz.squaredNorm() / m.squaredNorm();
        double l2 = 1 / (sin2_psi / pow(axis_l_.x(), 2) + cos2_psi / pow(axis_l_.z(), 2));
        double l = sqrt(l2);
        double d = (pos - p).norm();
        *dis = d - l;
        if (g) {
            double dis1, dis2, dis3;
            double del = 1e-4;
            get_dis_ellipsoid(pos + Vector3d(del, 0, 0), t, &dis1);
            get_dis_ellipsoid(pos + Vector3d(0, del, 0), t, &dis2);
            get_dis_ellipsoid(pos + Vector3d(0, 0, del), t, &dis3);
            *g = Vector3d((dis1 - *dis) / del, (dis2 - *dis) / del, (dis3 - *dis) / del);
            // cout << g->transpose() << endl;
        }
    }

    void get_dis_ellipsoid2(const Vector3d pos, const double t, double *dis, Vector3d *g = nullptr) const {
        double r = pos_ratio_ + dir_ * vel_ratio_ * t;
        Vector2d _p = r * start_pos_ + (1 - r) * end_pos_;
        Vector3d p(_p.x(), _p.y(), 1.5);
        Vector3d m = pos - p;
        *dis = pow(m.x(), 2) / pow(axis_l_.x() + 0.4, 2)
                + pow(m.y(), 2) / pow(axis_l_.y() + 0.4, 2)
                + pow(m.z(), 2) / pow(axis_l_.z() + 0.4, 2)
                - 1 + 0.4;
        if (g) {
            double dis1, dis2, dis3;
            double del = 1e-10;
            get_dis_ellipsoid2(pos + Vector3d(del, 0, 0), t, &dis1);
            get_dis_ellipsoid2(pos + Vector3d(0, del, 0), t, &dis2);
            get_dis_ellipsoid2(pos + Vector3d(0, 0, del), t, &dis3);
            *g = Vector3d((dis1 - *dis) / del, (dis2 - *dis) / del, (dis3 - *dis) / del);
            // cout << g->transpose() << endl;
        }
    }
};