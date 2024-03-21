#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>
#include <Eigen/Dense>
#include <chrono>

#include "bspline/uniform_bspline.hpp"
#include "map/sdf.hpp"

using namespace std;
using namespace Eigen;

class BsplineOptimizer2 {
public:
    BsplineOptimizer2();
    ~BsplineOptimizer2();
    int optimize(MatrixXd &ctrl_pts
        , const Vector3d &start_p, const Vector3d &start_v
        , const Vector3d &end_p, const Vector3d &end_v
        , const SdfMap *const sdf, const double ts
        , const double lamda_s, const double lamda_c, const double lamda_v
        , const double lamda_a, const double lamda_l, const double lamda_ep, const double lamda_ev
        , const double risk_dis, const double vmax
        , const double amax, vector<tuple<double, double, double, double, double, double>> &cost_history);
    
private:
    double lamda_s_, lamda_c_, lamda_v_, lamda_a_, lamda_l_, lamda_ep_, lamda_ev_, risk_dis_, vmax_, amax_, ts_;
    const SdfMap *sdf_map_;
    Vector3d start_p_, start_v_, end_p_, end_v_;
    chrono::steady_clock::time_point time_;
    int cnt_;
    vector<tuple<double, double, double, double, double, double>> cost_history_;
    static double cost_func(const vector<double> &u, vector<double> &grad, BsplineOptimizer2 *instance);
};