#pragma once

#include <vector>

#include "map/map.hpp"

class KinodynamicAstar {
public:
    class NodeInfo {
    public:
        double g_cost_;
        double f_cost_;
        Vector3d pos_;
        Vector3d vel_;
        Vector3d input_;
        double duration_;
        double time_;
        NodeInfo *parent_;
        char *visited_;
        Vector3i grid_idx_;
        NodeInfo() : g_cost_(0.0), f_cost_(0.0), time_(0.0), parent_(nullptr), visited_(0), grid_idx_(0, 0, 0) {};
    };
    typedef NodeInfo* NodeInfoPtr;
    class SearchNode {
    public:
        class SearchNodeComparator {
        public:
            bool operator()(SearchNode *node1, SearchNode *node2) {
                return node1->f_cost_ > node2->f_cost_;
            }
        };
    public:
        Vector3i grid_idx_;
        NodeInfo &info_;
        /*
        * 拷贝值用于排序（重要），如果使用info.f_cost_进行排序，
        * 则会存在bug（当修改openset中元素值时，priority_queue排序会存在问题）。
        */
        double g_cost_;
        double f_cost_;

        SearchNode(Vector3i idx, NodeInfo &info, double g_cost, double f_cost
            , Vector3d pos, Vector3d vel, Vector3d input
            , double duration, double time, NodeInfoPtr parent)
            : grid_idx_(idx), info_(info) {
            g_cost_ = g_cost;
            f_cost_ = f_cost;
            info_.g_cost_ = g_cost;
            info_.f_cost_ = f_cost;
            info_.pos_ = pos;
            info_.vel_ = vel;
            info_.input_ = input;
            info_.duration_ = duration;
            info_.time_ = time;
            info_.parent_ = parent;
        }
    };
    typedef SearchNode* SearchNodePtr;
private:
    GridMap &gridmap_;

    //参数
    Vector3d origin_, map_size_, resolution_;
    Vector3i map_grid_size_;
    double w_time_, max_vel_, max_acc_, tie_breaker_;
    double acc_resolution_, time_resolution_;
    double max_duration_;
    double safety_check_res_;
    double lambda_heu_;

    MatrixXd coef_shot_;
    double t_shot_;
    bool is_shot_succ_;
    vector<NodeInfo> path_;
    NodeInfo *nodeinfos;
    char *nodevisited_;

public:
    KinodynamicAstar(GridMap &map) : gridmap_(map) {
        origin_ = Vector3d(0., 0., 0.);
        resolution_ = Vector3d(1., 1., 1.) * gridmap_.resolution();
        map_grid_size_ = gridmap_.isize();
        map_size_ = map_grid_size_.cast<double>() * gridmap_.resolution();
        nodeinfos = new NodeInfo[map_grid_size_(2) * map_grid_size_(1) * map_grid_size_(0)];
        nodevisited_ = new char[map_grid_size_(2) * map_grid_size_(1) * map_grid_size_(0)];
#define nodeinfos(x, y, z) nodeinfos[z * map_grid_size_(1) * map_grid_size_(0) + y * map_grid_size_(0) + x]
#define nodevisited_(x, y, z) nodevisited_[z * map_grid_size_(1) * map_grid_size_(0) + y * map_grid_size_(0) + x]
        for (int i = 0; i < map_grid_size_(2); i++) {
            for (int j = 0; j < map_grid_size_(1); j++) {
                for (int k = 0; k < map_grid_size_(0); k++) {
                    nodeinfos(k, j, i).grid_idx_ = Vector3i(k, j, i);
                    nodeinfos(k, j, i).f_cost_ = 0;
                    nodeinfos(k, j, i).g_cost_ = 0;
                    nodeinfos(k, j, i).time_ = 0;
                    nodeinfos(k, j, i).parent_ = nullptr;
                    nodeinfos(k, j, i).visited_ = &nodevisited_(k, j, i);
                }
            }
        }
    }
    ~KinodynamicAstar() {
        delete nodeinfos;
        delete nodevisited_;
    }
    void set_param(
        double w_time, double max_vel, double max_acc, double tie_breaker
        , double acc_resolution, double time_resolution, double max_duration
        , double safety_check_res, double lambda_heu) {
        w_time_ = w_time;
        max_vel_ = max_vel;
        max_acc_ = max_acc;
        tie_breaker_ = tie_breaker;
        acc_resolution_ = acc_resolution;
        time_resolution_ = time_resolution;
        max_duration_ = max_duration;
        safety_check_res_ = safety_check_res;
        lambda_heu_ = lambda_heu;
    }
    bool search(Vector3d start_pt, Vector3d start_v,
        Vector3d end_pt, Vector3d end_v,
        const vector<DynObs> *dynobs=nullptr);
    pair<vector<Vector3d>, vector<Vector3d>> get_sample_path(double &dt);

private:
    vector<double> cubic(double a, double b, double c, double d);
    vector<double> quartic(double a, double b, double c, double d, double e);
    double estimateHeuristic(const Vector3d &p0, const Vector3d &v0
        , const Vector3d &p1, const Vector3d &v1, double& optimal_time);
    bool computeShotTraj(VectorXd state1, VectorXd state2, double time_to_goal);
    bool getGridOcc(const Vector3d &pos) {return gridmap_(pos(0), pos(1), pos(2));}
    Vector3d idx2pos(const Vector3i &idx) {
        return Vector3d((idx(0) + 0.5) * resolution_(0)
            , (idx(1) + 0.5) * resolution_(1), (idx(2) + 0.5) * resolution_(2));
    }
    Vector3i pos2idx(const Vector3d &pos) {
        Vector3i tmp = Vector3i(pos(0) / resolution_(0)
            , pos(1) / resolution_(1), fabs(resolution_(2)) > 1e-6 ? pos(2) / resolution_(2) : 0);
        if (tmp.x() < 0) {
            tmp.x() = 0;
        } else if (tmp.x() >= map_grid_size_.x()) {
            tmp.x() = map_grid_size_.x() - 1;
        }

        if (tmp.y() < 0) {
            tmp.y() = 0;
        } else if (tmp.y() >= map_grid_size_.y()) {
            tmp.y() = map_grid_size_.y() - 1;
        }

        if (tmp.z() < 0) {
            tmp.z() = 0;
        } else if (tmp.z() >= map_grid_size_.z()) {
            tmp.z() = map_grid_size_.z() - 1;
        }

        return tmp;
    }
    bool check_in_map(const Vector3d &pos) {
        if (pos(0) < origin_(0) || pos(1) < origin_(1) || pos(2) < origin_(2)
            || pos(0) > origin_(0) + map_size_(0) || pos(1) > origin_(1) + map_size_(1)
            || pos(2) > origin_(2) + map_size_(2)) {
            return false;
        }
        return true;
    }
    void stateTransit(const Vector3d &p0, const Vector3d &v0
        , Vector3d &p1, Vector3d &v1, Vector3d &um, const double &tau)
    {
        p1 = p0 + 0.5 * tau * tau * um + v0 * tau ;
        v1 = v0 + tau * um;
    }
};