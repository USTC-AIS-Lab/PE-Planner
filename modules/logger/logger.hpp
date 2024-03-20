#pragma once

#include <iostream>
#include <fstream>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

class Logger {
public:
    ofstream file_;

    //data
    Vector3d pos_;
    Vector3d vel_;
    Vector3d acc_;
    Vector3d att_;
    Vector3d rate_;
    double tilt_;
    Vector3d disturbance_;
    Vector4d u_;
    double dis_to_obs_;
    double solution_time_;
    double length_;
    double time_;

    Logger(string logname) {
        time_t timep;
        time(&timep);
        char tmp[64];
        strftime(tmp, sizeof(tmp), ("%Y-%m-%d-%H-%M-%S--" + logname + ".txt").c_str(), localtime(&timep));
        file_.open(tmp, fstream::out);
        if (!file_) {
            cerr << "Fail to open map file" << endl;
            throw "Fail to open logger file";
        }
        cout << "Create log file " << tmp << endl;
    }
    ~Logger() {
        file_.close();
    }

    void update() {
        file_ << pos_(0) << " " << pos_(1) << " " << pos_(2) << " ";
        file_ << vel_(0) << " " << vel_(1) << " " << vel_(2) << " ";
        file_ << acc_(0) << " " << acc_(1) << " " << acc_(2) << " ";
        file_ << att_(0) << " " << att_(1) << " " << att_(2) << " ";
        file_ << rate_(0) << " " << rate_(1) << " " << rate_(2) << " ";
        file_ << tilt_ << " ";
        file_ << disturbance_(0) << " " << disturbance_(1) << " " << disturbance_(2) << " ";
        file_ << u_(0) << " " << u_(1) << " " << u_(2) << " " << u_(3) << " ";
        file_ << dis_to_obs_ << " ";
        file_ << solution_time_ << " ";
        file_ << length_ << " ";
        file_ << time_ << " ";
        file_ << endl;
    }
};