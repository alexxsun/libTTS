// base code includes
#include "../ANN/ANN/ANN.h"

#include <map>
#include <string>
#include <vector>

using namespace std;

int als_segment(const string &infile, double th_radius, double th_forest_pers, double th_height, bool th_height_funtype,
            double th_pers_H, double th_pers_I, const string &ds);

int read_points(const string &infile, vector<vector<double>> &pts);

int connect_projected_pts(const vector<vector<double>> &pts, double &Radius, vector<double> &density,
                          vector<vector<int>> &VV);

int segment_forest(const double &perc_pers, const vector<double> &density, vector<vector<int>> &VV,
                   vector<int> &pts_lbls, map<int, vector<int>> &tree_clusters);

int compute_cylinder(double &th_height, const double &radius, const vector<vector<double>> &pts,
                     const vector<int> &pids, vector<vector<int>> &VV, double &zmin, double &zmax);

int segment_cluster(const vector<vector<double>> &pts, const double &zmin, const double &zmax,
                    const vector<int> &vertices, const vector<vector<int>> &VV, const bool &height_funtype,
                    const double &perc_pers, vector<int> &pts_lbls, int &regionN);

int merge_labels(const vector<int> &lbls1, const vector<int> &lbls2, vector<int> &lbls, int &trees_num);

int post_process(vector<int> &pts_lbls, const int &trees_num, const string &ds);

void output_labeled_pts(vector<vector<double>> &points, vector<int> &pids, const string &outfile);
