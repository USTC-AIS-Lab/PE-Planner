#pragma once

#include <iostream>
#include <limits>
#include <math.h>
#include <Eigen/Dense>
#include <chrono>

#include "map.hpp"

using namespace std;
using namespace Eigen;

class SdfMap
{
public:
    double resolution_;
    Vector3d map_size_;
    Vector3i voxel_num_;
    vector<double> sdf_map_;
public:
    SdfMap(double resolution, GridMap &grid_map);
    ~SdfMap();
    inline int to_address(int x, int y, int z) const {return z * voxel_num_(1) * voxel_num_(0) + y * voxel_num_(0) + x;}
    inline double get_dis(Vector3i &idx) const {
        if (check_in_map(idx))
            return sdf_map_[idx.z() * voxel_num_(1) * voxel_num_(0) + idx.y() * voxel_num_(0) + idx.x()];
        else
            return -5.0;
    }
    template <typename T> inline pair<T, Matrix<T, 3, 1>> get_dist_with_grad_trilinear(Matrix<T, 3, 1> pos) const;
    bool check_in_map(Vector3d &pos) const {
        if (pos.x() > map_size_(0) || pos.y() > map_size_(1) || pos.z() > map_size_(2) 
            || pos.x() < 0 || pos.y() < 0 || pos.z() < 0) {
            return false;
        }
        return true;
    }
    bool check_in_map(Vector3i &idx) const {
        if (idx.x() >= voxel_num_(0) || idx.y() >= voxel_num_(1) || idx.z() >= voxel_num_(2) 
            || idx.x() < 0 || idx.y() < 0 || idx.z() < 0) {
            return false;
        }
        return true;
    }
    const Vector3i &voxel_num() const {return voxel_num_;}
private:
    inline int to_address_(int x, int y, int z) const {return z * voxel_num_(1) * voxel_num_(0) + y * voxel_num_(0) + x;}
    template <typename F_get_val, typename F_set_val>
    inline void fill_esdf(F_get_val f_get_val, F_set_val f_set_val, int start, int end) {
        int v[end - start + 1];
        double z[end - start + 2];

        int k = start;
        v[start] = start;
        z[start] = -std::numeric_limits<double>::max();
        z[start + 1] = std::numeric_limits<double>::max();

        for (int q = start + 1; q <= end; q++) {
            k++;
            double s;

            do {
            k--;
            s = ((f_get_val(q) + q * q) - (f_get_val(v[k]) + v[k] * v[k])) / (2 * q - 2 * v[k]);
            } while (s <= z[k]);

            k++;

            v[k] = q;
            z[k] = s;
            z[k + 1] = std::numeric_limits<double>::max();
        }

        k = start;

        for (int q = start; q <= end; q++) {
            while (z[k + 1] < q) k++;
            double val = (q - v[k]) * (q - v[k]) + f_get_val(v[k]);
            f_set_val(q, val);
        }
    }
};

inline SdfMap::SdfMap(double resolution, GridMap &grid_map)
{
    auto t0 = chrono::steady_clock::now();
    resolution_ = resolution;
    map_size_ = grid_map.size();
    voxel_num_ = Vector3i(int(grid_map.size()(0) / resolution_)
        , int(grid_map.size()(1) / resolution_)
        , int(grid_map.size()(2) / resolution_));
    vector<double> tmp_buf1(voxel_num_(0) * voxel_num_(1) * voxel_num_(2));
    vector<double> tmp_buf2(voxel_num_(0) * voxel_num_(1) * voxel_num_(2));
    vector<double> pdis_buf(voxel_num_(0) * voxel_num_(1) * voxel_num_(2));
    vector<double> ndis_buf(voxel_num_(0) * voxel_num_(1) * voxel_num_(2));
    double resolution_ratio = resolution_ / grid_map.resolution();

    //positive
    for (int x = 0; x < voxel_num_(0); x++) {
        for (int y = 0; y < voxel_num_(1); y++) {
            fill_esdf(
                [&](int z) {
                    return grid_map(int(x * resolution_ratio), int(y * resolution_ratio), int(z * resolution_ratio)) == 1 ?
                        0.0 :
                        numeric_limits<double>::max();
                }
                , [&](int z, double val) {tmp_buf1[to_address_(x, y, z)] = val;}
                , 0, voxel_num_(2) - 1);
        }
    }

    for (int x = 0; x < voxel_num_(0); x++) {
        for (int z = 0; z < voxel_num_(2); z++) {
            fill_esdf(
                [&](int y) {return tmp_buf1[to_address_(x, y, z)];}
                , [&](int y, double val) {tmp_buf2[to_address_(x, y, z)] = val;}
                , 0, voxel_num_(1) - 1);
        }
    }

    for (int y = 0; y < voxel_num_(1); y++) {
        for (int z = 0; z < voxel_num_(2); z++) {
            fill_esdf(
                [&](int x) {return tmp_buf2[to_address_(x, y, z)];}
                , [&](int x, double val) {pdis_buf[to_address_(x, y, z)] = resolution_ * sqrt(val);}
                , 0, voxel_num_(0) - 1);
        }
    }

    //negative
    for (int x = 0; x < voxel_num_(0); x++) {
        for (int y = 0; y < voxel_num_(1); y++) {
            fill_esdf(
                [&](int z) {
                    return grid_map(int(x * resolution_ratio), int(y * resolution_ratio), int(z * resolution_ratio)) == 0 ?
                        0.0 :
                        numeric_limits<double>::max();
                }
                , [&](int z, double val) {tmp_buf1[to_address_(x, y, z)] = val;}
                , 0, voxel_num_(2) - 1);
        }
    }

    for (int x = 0; x < voxel_num_(0); x++) {
        for (int z = 0; z < voxel_num_(2); z++) {
            fill_esdf(
                [&](int y) {return tmp_buf1[to_address_(x, y, z)];}
                , [&](int y, double val) {tmp_buf2[to_address_(x, y, z)] = val;}
                , 0, voxel_num_(1) - 1);
        }
    }

    for (int y = 0; y < voxel_num_(1); y++) {
        for (int z = 0; z < voxel_num_(2); z++) {
            fill_esdf(
                [&](int x) {return tmp_buf2[to_address_(x, y, z)];}
                , [&](int x, double val) {ndis_buf[to_address_(x, y, z)] = resolution_ * sqrt(val);}
                , 0, voxel_num_(0) - 1);
        }
    }

    for (int z = 0; z < voxel_num_(2); z++) {
        for (int y = 0; y < voxel_num_(1); y++) {
            for (int x = 0; x < voxel_num_(0); x++) {
                if (ndis_buf[to_address_(x, y, z)] > 0.0)
                    sdf_map_.push_back(pdis_buf[to_address_(x, y, z)] - ndis_buf[to_address_(x, y, z)] + resolution_ / 2.0);
                else
                    sdf_map_.push_back(pdis_buf[to_address_(x, y, z)] - resolution_ / 2.0);
            }
        }
    }
    cerr << "Build ESDF, spend " << chrono::duration<double>(chrono::steady_clock::now() - t0).count() * 1e3 << "ms" << endl;
}

inline SdfMap::~SdfMap()
{

}

template <typename T>
inline pair<T, Matrix<T, 3, 1>> SdfMap::get_dist_with_grad_trilinear(Matrix<T, 3, 1> pos) const {
    Vector3i idx(int((pos.x() + 0.5 * resolution_) / resolution_) - 1
        , int((pos.y() + 0.5 * resolution_) / resolution_) - 1
        , int((pos.z() + 0.5 * resolution_) / resolution_) - 1);
    Matrix<T, 3, 1> idx_pos((idx.x() + 0.5) * resolution_
        , (idx.y() + 0.5) * resolution_
        , (idx.z() + 0.5) * resolution_);
    Matrix<T, 3, 1> diff = (pos - idx_pos) / resolution_;
    T values[2][2][2];
    for (int x = 0; x < 2; x++) {
        for (int y = 0; y < 2; y++) {
            for (int z = 0; z < 2; z++) {
                Vector3i current_idx = idx + Vector3i(x, y, z);
                values[x][y][z] = get_dis(current_idx);
            }
        }
    }
    T v00 = (1 - diff[0]) * values[0][0][0] + diff[0] * values[1][0][0];
    T v01 = (1 - diff[0]) * values[0][0][1] + diff[0] * values[1][0][1];
    T v10 = (1 - diff[0]) * values[0][1][0] + diff[0] * values[1][1][0];
    T v11 = (1 - diff[0]) * values[0][1][1] + diff[0] * values[1][1][1];
    T v0 = (1 - diff[1]) * v00 + diff[1] * v10;
    T v1 = (1 - diff[1]) * v01 + diff[1] * v11;
    T dis = (1 - diff[2]) * v0 + diff[2] * v1;

    Matrix<T, 3, 1> grad;
    grad[2] = (v1 - v0) / resolution_;
    grad[1] = ((1 - diff[2]) * (v10 - v00) + diff[2] * (v11 - v01)) / resolution_;
    grad[0] = (1 - diff[2]) * (1 - diff[1]) * (values[1][0][0] - values[0][0][0]);
    grad[0] += (1 - diff[2]) * diff[1] * (values[1][1][0] - values[0][1][0]);
    grad[0] += diff[2] * (1 - diff[1]) * (values[1][0][1] - values[0][0][1]);
    grad[0] += diff[2] * diff[1] * (values[1][1][1] - values[0][1][1]);
    grad[0] /= resolution_;
    return make_pair(dis, grad);
}
