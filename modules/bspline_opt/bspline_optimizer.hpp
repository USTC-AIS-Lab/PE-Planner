#pragma once

#include <math.h>

#include <Eigen/Dense>
#include <chrono>

#include "map/sdf.hpp"

using namespace Eigen;

class BsplineOptimizer {
public:
    BsplineOptimizer() {}
    ~BsplineOptimizer() {}
    int optimize(MatrixXd &ctrl_pts
        , const Vector3d &start_p, const Vector3d &start_v, const Vector3d &start_a
        , const Vector3d &end_p, const Vector3d &end_v, const Vector3d &end_a
        , const SdfMap *const sdf, const double ts
        , const double lamda_s, const double lamda_c, const double lamda_v
        , const double lamda_a, const double lamda_l, const double lamda_dl
        , const double lamda_ep, const double lamda_ev, const double lamda_ea
        , const double risk_dis, const double vmax
        , const double amax, vector<tuple<double, double, double, double, double, double>> &cost_history
        , const vector<DynObs> *dynobs=nullptr);
        
private:
    double lamda_s_, lamda_c_, lamda_v_, lamda_a_, lamda_l_, lamda_dl_, lamda_ep_, lamda_ev_, lamda_ea_, risk_dis_, vmax_, amax_, ts_;
    const SdfMap *sdf_map_;
    const vector<DynObs> *dynobs_;
    MatrixXd *ctrl_pts_;
    Vector3d start_p_, start_v_, start_a_, end_p_, end_v_, end_a_;
    chrono::steady_clock::time_point time_;
    int cnt_;
    vector<tuple<double, double, double, double, double, double>> cost_history_;
    static double cost_function(BsplineOptimizer *instance, const VectorXd &x, VectorXd &g);
};
