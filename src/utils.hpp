#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "map/map.hpp"

using namespace std;
using namespace Eigen;

inline double get_dis_to_dynobs(const vector<DynObs> &dynobs, const Vector3d &pos) {
    double dis = INFINITY;
    for (const auto &o : dynobs) {
        double tmp;
        o.get_dis(pos, 0, &tmp);
        dis = min(dis, tmp);
    }
    return dis;
}