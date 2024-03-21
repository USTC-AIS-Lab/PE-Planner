#include <iostream>
#include <vector>
#include <queue>
#include <Eigen/Dense>
#include <chrono>
#include <iomanip>

#include "kinodynamic_astar.hpp"

vector<double> KinodynamicAstar::cubic(double a, double b, double c, double d) {
    vector<double> dts;

    double a2 = b / a;
    double a1 = c / a;
    double a0 = d / a;

    double Q = (3 * a1 - a2 * a2) / 9;
    double R = (9 * a1 * a2 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
    double D = Q * Q * Q + R * R;
    if (D > 0) {
        double S = std::cbrt(R + sqrt(D));
        double T = std::cbrt(R - sqrt(D));
        dts.push_back(-a2 / 3 + (S + T));
        return dts;
    } else if (D == 0) {
        double S = std::cbrt(R);
        dts.push_back(-a2 / 3 + S + S);
        dts.push_back(-a2 / 3 - S);
        return dts;
    } else {
        double theta = acos(R / sqrt(-Q * Q * Q));
        dts.push_back(2 * sqrt(-Q) * cos(theta / 3) - a2 / 3);
        dts.push_back(2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) - a2 / 3);
        dts.push_back(2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) - a2 / 3);
        return dts;
    }
}

vector<double> KinodynamicAstar::quartic(double a, double b, double c, double d, double e) {
    vector<double> dts;

    double a3 = b / a;
    double a2 = c / a;
    double a1 = d / a;
    double a0 = e / a;

    vector<double> ys = cubic(1, -a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);
    double y1 = ys.front();
    double r = a3 * a3 / 4 - a2 + y1;
    if (r < 0)
        return dts;

    double R = sqrt(r);
    double D, E;
    if (R != 0) {
        D = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 + 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
        E = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 - 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
    } else {
        D = sqrt(0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0));
        E = sqrt(0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0));
    }

    if (!std::isnan(D)) {
        dts.push_back(-a3 / 4 + R / 2 + D / 2);
        dts.push_back(-a3 / 4 + R / 2 - D / 2);
    }
    if (!std::isnan(E)) {
        dts.push_back(-a3 / 4 - R / 2 + E / 2);
        dts.push_back(-a3 / 4 - R / 2 - E / 2);
    }

    return dts;
}

double KinodynamicAstar::estimateHeuristic(const Vector3d &p0, const Vector3d &v0
    , const Vector3d &p1, const Vector3d &v1, double& optimal_time) {
    const Vector3d dp = p1 - p0;

    double c1 = -36 * dp.dot(dp);
    double c2 = 24 * (v0 + v1).dot(dp);
    double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
    double c4 = 0;
    double c5 = w_time_;

    std::vector<double> ts = quartic(c5, c4, c3, c2, c1);

    double v_max = max_vel_ * 0.5;
    double t_bar = (p0 - p1).lpNorm<Infinity>() / v_max;
    ts.push_back(t_bar);

    double cost = 100000000;
    double t_d = t_bar;

    for (auto t : ts) {
        if (t < t_bar)
            continue;
        double c = -c1 / (3 * t * t * t) - c2 / (2 * t * t) - c3 / t + w_time_ * t;
        if (c < cost) {
            cost = c;
            t_d = t;
        }
    }

    optimal_time = t_d;

    return 1.0 * (1 + tie_breaker_) * cost;
}

bool KinodynamicAstar::computeShotTraj(Eigen::VectorXd state1, Eigen::VectorXd state2, double time_to_goal) {
    /* ---------- get coefficient ---------- */
    const Vector3d p0 = state1.head(3);
    const Vector3d dp = state2.head(3) - p0;
    const Vector3d v0 = state1.segment(3, 3);
    const Vector3d v1 = state2.segment(3, 3);
    const Vector3d dv = v1 - v0;
    double t_d = time_to_goal;
    MatrixXd coef(3, 4);

    Vector3d a = 1.0 / 6.0 * (-12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv);
    Vector3d b = 0.5 * (6.0 / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv);
    Vector3d c = v0;
    Vector3d d = p0;

    // 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0
    // a*t^3 + b*t^2 + v0*t + p0
    coef.col(3) = a, coef.col(2) = b, coef.col(1) = c, coef.col(0) = d;

    Vector3d coord, vel, acc;
    VectorXd poly1d, t, polyv, polya;
    Vector3i index;

    Eigen::MatrixXd Tm(4, 4);
    Tm << 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0;

    /* ---------- forward checking of trajectory ---------- */
    double t_delta = t_d / 20;
    for (double time = t_delta; time <= t_d; time += t_delta) {
        t = VectorXd::Zero(4);
        for (int j = 0; j < 4; j++)
            t(j) = pow(time, j);

        for (int dim = 0; dim < 3; dim++) {
            poly1d = coef.row(dim);
            coord(dim) = poly1d.dot(t);
            vel(dim) = (Tm * poly1d).dot(t);
            acc(dim) = (Tm * Tm * poly1d).dot(t);
        }

        if (coord(0) < origin_(0) || coord(0) >= map_size_(0) || coord(1) < origin_(1) || coord(1) >= map_size_(1) ||
            coord(2) < origin_(2) || coord(2) >= map_size_(2)) {
            return false;
        }

        if (getGridOcc(coord) == 1) {
            return false;
        }
    }
    coef_shot_ = coef;
    t_shot_ = t_d;
    is_shot_succ_ = true;
    return true;
}

bool KinodynamicAstar::search(Vector3d start_pt, Vector3d start_v,
    Vector3d end_pt, Vector3d end_v,
    const vector<DynObs> *dynobs) {
    auto t0 = chrono::steady_clock::now();
    memset(nodevisited_, 0, map_grid_size_(2) * map_grid_size_(1) * map_grid_size_(0));
    priority_queue<SearchNodePtr, vector<SearchNodePtr>, SearchNode::SearchNodeComparator> openset;
    vector<SearchNodePtr> nodeset;
    Vector3i end_idx = pos2idx(end_pt);
    Vector3i start_idx = pos2idx(start_pt);

    SearchNodePtr node = new SearchNode(
        start_idx, nodeinfos(start_idx.x(), start_idx.y(), start_idx.z())
        , 0, 0, start_pt, start_v, Vector3d(0., 0., 0.), 0, 0, nullptr);
    *nodeinfos(start_idx.x(), start_idx.y(), start_idx.z()).visited_ = 1;
    nodeset.push_back(node);
    openset.push(node);

    SearchNodePtr end_node = nullptr;

    int cnt = 0;
    
    while (!openset.empty()) {
        SearchNodePtr node;
        do {
            if (openset.empty()) {
                goto FAILED;
            }
            node = openset.top();
            openset.pop();
        } while(*node->info_.visited_ != 1);
        *node->info_.visited_ = -1;
        double time = node->info_.time_;
        Vector3d deltap = end_pt - node->info_.pos_;
        Vector3d deltav = end_v - node->info_.vel_;

        bool is_nearest = false;
        
        if (deltap.squaredNorm() < pow(0.6, 2)
            /*fabs(deltap(0)) < 0.2 && fabs(deltap(1)) < 0.2 && fabs(deltap(2)) < 0.2*/
            /*&& fabs(deltav(0)) < 1.0 && fabs(deltav(1)) < 1.0 && fabs(deltav(2)) < 1.0*/) {
            // double time_to_goal = 0.;
            // Matrix<double, 6, 1> state1, state2;
            // state1 << node->info_.pos_, node->info_.vel_;
            // state2 << end_pt, end_v;
            // estimateHeuristic(node->info_.pos_, node->info_.vel_, end_pt, end_v, time_to_goal);
            // is_shot_succ_ = false;
            // if (computeShotTraj(state1, state2, time_to_goal))
            //     is_nearest = true;
            is_nearest = true;
        }

        if (is_nearest) {
            auto t1 = chrono::steady_clock::now();
            // cerr << "Kinodynamic A* find avaliable path, spend " << chrono::duration<double>(t1 - t0).count() * 1e3 << "ms" << ", final cost: " << node->f_cost_ << endl;
            end_node = node;
            break;
        }

        //运动原语
        vector<Vector3d> inputs;
        vector<double> durations;
        for (double ax = -max_acc_; ax < max_acc_ + 1e-7; ax += max_acc_ * acc_resolution_) {
            for (double ay = -max_acc_; ay < max_acc_ + 1e-7; ay += max_acc_ * acc_resolution_) {
                for (double az = -max_acc_; az < max_acc_ + 1e-7; az += max_acc_ * acc_resolution_) {
                    inputs.push_back(Vector3d(ax, ay, az));
                }
            }
        }
        for (double t = max_duration_ * time_resolution_; t < max_duration_ + 1e-6; t += max_duration_ * time_resolution_) {
            durations.push_back(t);
        }

        vector<NodeInfoPtr> tmp_expand_nodes;

        for (auto &input : inputs) {
            for (auto &tau : durations) {
                //状态转移
                Vector3d new_pos, new_vel;
                stateTransit(node->info_.pos_, node->info_.vel_, new_pos, new_vel, input, tau);
                double new_time = time + tau;

                if (!check_in_map(new_pos)) {
                    continue;
                }
                
                //检查是否在closeset
                Vector3i new_idx = pos2idx(new_pos);
                if (new_idx == node->grid_idx_) {
                    continue;
                }
                if (*nodeinfos(new_idx(0), new_idx(1), new_idx(2)).visited_ == -1) {
                    continue;
                }

                //安全检查
                if ((fabs(new_vel(0)) > max_vel_ || fabs(new_vel(1)) > max_vel_ || fabs(new_vel(2)) > max_vel_)) {
                    continue;
                }
                bool safe = true;
                for (double t = tau * safety_check_res_; t < tau + 1e-6; t += tau * safety_check_res_) {
                    Vector3d p, v;
                    stateTransit(node->info_.pos_, node->info_.vel_, p, v, input, t);
                    if (getGridOcc(p)) {
                        safe = false;
                        break;
                    }
                    if (dynobs) {
                        for (auto &o : *dynobs) {
                            double dis = 0.0;
                            o.get_dis_ellipsoid(p, time + t, &dis);
                            if (dis < 0) {
                                safe = false;
                                break;
                            }
                        }
                    }
                }
                if (!safe) {
                    continue;
                }
                
                //计算代价
                double gcost = node->g_cost_ + (input.squaredNorm() + w_time_) * tau;
                double optimal_time = 0;
                double fcost = gcost + lambda_heu_ * estimateHeuristic(new_pos, new_vel, end_pt, end_v, optimal_time);

                bool prune = false;
                for (auto &expand_node : tmp_expand_nodes) {
                    if (expand_node->grid_idx_ == new_idx) {
                        prune = true;
                        if (fcost < expand_node->f_cost_) {
                            expand_node->f_cost_ = fcost;
                            expand_node->g_cost_ = gcost;
                            expand_node->input_ = input;
                            expand_node->pos_ = new_pos;
                            expand_node->vel_ = new_vel;
                            expand_node->duration_ = tau;
                            expand_node->time_ = new_time;
                        }
                        break;
                    }
                }

                if (prune) {
                    continue;
                }

                if (*nodeinfos(new_idx(0), new_idx(1), new_idx(2)).visited_ == 0) {
                    SearchNodePtr new_node = new SearchNode(new_idx
                        , nodeinfos(new_idx(0), new_idx(1), new_idx(2))
                        , gcost, fcost, new_pos, new_vel, input, tau, new_time
                        , &(node->info_));
                    *new_node->info_.visited_ = 1;
                    nodeset.push_back(new_node);
                    openset.push(new_node);
                    tmp_expand_nodes.push_back(&(new_node->info_));
                } else if (*nodeinfos(new_idx(0), new_idx(1), new_idx(2)).visited_ == 1) {
                    if (gcost < nodeinfos(new_idx(0), new_idx(1), new_idx(2)).g_cost_) {
                        SearchNodePtr new_node = new SearchNode(new_idx
                        , nodeinfos(new_idx(0), new_idx(1), new_idx(2))
                        , gcost, fcost, new_pos, new_vel, input, tau, new_time
                        , &(node->info_));
                        nodeset.push_back(new_node);
                        openset.push(new_node);
                    }
                }
            }
        }

        if (++cnt > 400000) {
            cout << "Spend too much time" << endl;
            end_node = nullptr;
            break;
        }
    }
    if (end_node == nullptr) {
FAILED:
        auto t1 = chrono::steady_clock::now();
        for (int i = 0; i < nodeset.size(); i++) {
            delete nodeset[i];
        }
        // delete []nodeinfos;
        cerr << "Kinodynamic A* failed, spend " << chrono::duration<double>(t1 - t0).count() * 1e3 << "ms" << endl;
        return false;
    }
    NodeInfoPtr start_node = &(end_node->info_);
    path_.clear();
    if (start_node != nullptr)
        path_.push_back(*start_node);
    while (start_node != nullptr && start_node->parent_ != nullptr) {
        start_node = start_node->parent_;
        path_.push_back(*start_node);
    }
    reverse(path_.begin(), path_.end());
    for (int i = 0; i < nodeset.size(); i++) {
        delete nodeset[i];
    }
    // delete []nodeinfos;

    return true;
}

pair<vector<Vector3d>, vector<Vector3d>> KinodynamicAstar::get_sample_path(double &dt) {
    vector<Vector3d> sample_path, sample_vel;
    double total_t = path_[path_.size() - 1].time_;
    double num = ceil(total_t / dt);
    dt = total_t / num;
    double t = 0.0;
    int idx = 0;
    do {
        Vector3d p0 = path_[idx].pos_;
        Vector3d v0 = path_[idx].vel_;
        Vector3d u = path_[idx + 1].input_;
        Vector3d p1, v1;
        stateTransit(p0, v0, p1, v1, u, t - path_[idx].time_);
        sample_path.push_back(p1);
        sample_vel.push_back(v1);
        t += dt;
        while (idx + 1 < path_.size() - 1 && t >= path_[idx + 1].time_) {
            idx++;
        }
    } while (t < path_[path_.size() - 1].time_ + 1e-6);
    for (t = 0; t < t_shot_; t += dt) {
        Vector3d pos, vel;
        Vector4d poly1d, time, time2;

        for (int j = 0; j < 4; j++) {
            time(j) = pow(t, j);
            if (j >= 1)
                time2(j) = j * pow(t, j - 1);
            else
                time2(j) = 0;
        }

        for (int dim = 0; dim < 3; dim++)
        {
            poly1d = coef_shot_.row(dim);
            pos(dim) = poly1d.dot(time);
            vel(dim) = poly1d.dot(time2);
        }
        sample_path.push_back(pos);
        sample_vel.push_back(vel);
    }

    return make_pair(sample_path, sample_vel);
}