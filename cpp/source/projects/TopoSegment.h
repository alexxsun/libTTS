//
// Created by alex on 12/26/23.
// Core class for Topology-based Tree Segmentation.
//

#ifndef TTS_TOPOSEGMENT_H
#define TTS_TOPOSEGMENT_H

#include "../forman/formangradient.h"
#include <filesystem> // C++ 17

using namespace std;

class TopoSegment : public FormanGradient {
protected:
    bool _showlog = false;

    std::map<implicitS, int> _min2lbl;
    // min -> label. valid label starts from 1, <0: invalid labels. 0: not labeled.

    std::vector<int> _pts_lbls;
    size_t _segs_num = 0;

    int _init() {
        if (sc.getVertexCoordSize() != 4) {
            return 1; // sth wrong
        }
        // compute the Forman gradient
        computeFormanGradient(true);
        return 0; // everything goes well
    }

    // generate vv relationship of points (x,y,z,...) with p.z in [th_h1, th_h2]
    void _get_connected_points(const double &th_h1, const double &th_h2, std::vector<std::set<int> > &con_pts);

    // label all points from critical mins: -2: boundary, 0: unlabeled, >=1: valid labels
    std::vector<int> _label_points_by_grouped_mins(const std::map<int, std::vector<implicitS> > &gp_mins);

    // Todo: parallel version of _label_points_by_grouped_mins
    int _label_points_by_grouped_mins_parallel(const std::map<int, std::vector<implicitS> > &gp_mins);

    // write critical mins: xyz
    void _output_mins_pts(const std::vector<implicitS> &mins, const string &ptsfile, const bool &scaled);

    // write critical mins: xyz in vtk format
    void _output_mins_vtk(const std::vector<implicitS> &mins, const string &vtktfile, const bool &scaled);

    // write grouped critical mins (x,y,z,gid)
    void _output_gpmins_pts(const map<int, std::vector<implicitS> > &gp_mins, const string &ptsfile,
                            const bool &scaled);

    // write labeled points. xyzl.
    void _output_pts_with_label_pts(const string &outfile, const std::vector<int> &lbls, const bool &scaled);

    void _output_pts_with_label_pts_ply(const string &outfile, const std::vector<int> &lbls, const bool &scaled);

    // write labeled points. xyzl. in vtk format
    void _output_pts_with_label_vtk(const string &outfile, const std::vector<int> &lbls, const bool &scaled);

public:
    // constructor:
    TopoSegment(const string &infile, const int &fid, const bool &in_debug = false) : FormanGradient(infile, fid) {
        _showlog = in_debug;
        if (_init() > 0) {
            cerr << "Check infile!\nformat: xyzv\n";
            exit(1);
        }
    }

    // point clustering by Forman over-segmentation
    int cluster(string &out_ptsfile) {
        std::vector<implicitS> all_mins(criticalS[0].begin(), criticalS[0].end()); // all mins.
        std::map<implicitS, int> min2lbl;
        std::map<implicitS, std::set<implicitS> > min2mins; // critical minima to its directed connected critical minima

        // measure the time usage
        IO_Timer time;
        cout << "output over-segmentation results\n";
        std::map<int, std::vector<implicitS> > gp_mins;
        int lbl = 1; // min index of all_min = lbl-1
        for (const auto &m: all_mins) {
            gp_mins[lbl].push_back(m);
            min2lbl[m] = lbl;
            lbl = lbl + 1;
        }

        time.start();
        _label_points_by_grouped_mins_parallel(gp_mins);
        time.stop();
        cout << "pts labeled by mins computed " << time.getElapsedTime() << " s" << endl;

        if (out_ptsfile.empty()) {
            time.start();
            if (_file_extension == ".ply") {
                out_ptsfile = _workspace + _file_name + "_lbl.ply";
                _output_pts_with_label_pts_ply(out_ptsfile, _pts_lbls, false);
            } else {
                out_ptsfile = _workspace + _file_name + "_lbl.pts";
                _output_pts_with_label_pts(out_ptsfile, _pts_lbls, false);
            }
            time.stop();
            cout << "out_ptsfile=" << out_ptsfile << endl;
            cout << "pts labeled by mins written " << time.getElapsedTime() << " s" << endl;
        } else {
            time.stop();
        }

        return 0;
    }


    unsigned int get_points_num() {
        return sc.getVerticesNum();
    }

    std::vector<int> get_pts_label() {
        return _pts_lbls;
    }

    int get_segs_number() {
        // Constant member functions are those functions that are denied permission to
        // change the values of the data members of their class.
        // valid labels start from 1.
        _segs_num = 0;
        for (const auto &m: _min2lbl) {
            int lbl = m.second;
            if (lbl > 0) {
                _segs_num = _segs_num + 1;
            }
        }
        return _segs_num;
    }

    // initial min seeds are components of low-heights from input alpha shape
    // not suggested. keep here for the original paper.
    int initSeeds_mins_by_as();

    // init min labels from locations in 3D space.
    // location: tree bottom points, or tree trunk points.
    int initSeeds_mins_by_trunkpts_xyz(const std::vector<std::vector<float> > &seeds_pts,
                                       const double &th_p2trunk_distance = 0.2,
                                       const double &th_search_radius = 0.25
    );


    // tree growth: label minimums from initialized seeds
    int growFromSeeds_mins_basic();

    // label points from labeled mins via min's influence region
    string label_pts_from_mins(const bool &output = true);

    string label_pts_from_mins_parallel(const bool &output = true);

    // Output functions
    int write_results_xyzl_txt(std::string &outfile) {
        if (outfile.empty()) {
            outfile = _workspace + _file_name + "_lbl.pts";
        }
        cout << "write: " << outfile << endl;
        _output_pts_with_label_pts(outfile, _pts_lbls, false);

        return 0;
    }

    int write_results_xyzl_vtk(string &outfile) {
        if (outfile.empty()) {
            outfile = _workspace + _file_name + "_lbl.vtk";
        }
        cout << "write: " << outfile << endl;
        _output_pts_with_label_vtk(outfile, _pts_lbls, false);

        return 0;
    }

    string output_labeled_pts(const string &outfile);


    // Other functions
    int tts_read_pts_fs(const string &infile, std::vector<std::vector<float> > &ret_pts) {
        std::string line;
        ifstream fStream(infile);
        if (fStream.is_open()) {
            while (getline(fStream, line)) {
                if (line.empty()) {
                    continue;
                }
                istringstream iss(line);
                vector<float> coordinates;
                // https://stackoverflow.com/questions/9986091/how-do-you-convert-a-string-into-an-array-of-floats
                std::copy(std::istream_iterator<float>(iss), std::istream_iterator<float>(),
                          std::back_inserter(coordinates));
                ret_pts.push_back(coordinates);
            }
            fStream.close();
        }

        // debug
        std::cout << "pts#: " << ret_pts.size() << endl;
        if (!ret_pts.empty()) {
            std::cout << ret_pts[0][0] << " " << ret_pts[0][1] << " " << ret_pts[0][2] << " " << ret_pts[0][3] << "\n";
        }
        return 0; // everything goes well
    }

    float dis_sq_pt2pts(const std::array<float, 3> &pt, const std::set<std::array<float, 3> > &pts) {
        float dis = std::numeric_limits<float>::max();

        for (const auto &p: pts) {
            float dx, dy, dz;
            dx = pt[0] - p[0];
            dy = pt[1] - p[1];
            dz = pt[2] - p[2];

            float tmp_dis = dx * dx + dy * dy + dz * dz;
            if (tmp_dis < dis) {
                dis = tmp_dis;
            }
        }
        //dis = sqrt(dis);

        return dis;
    }
};


#endif //TTS_TOPOSEGMENT_H
