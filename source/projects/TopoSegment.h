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
    // Note:
    // This TopoSegment class doesn't do any postprocessing, e.g., validate and refine segmentation
    // to decouple steps/functions
    // TopoSegment just does the segmentation, doesn't care the final segmentation is reasonable or not
    // due to the data quality or user interests/settings of target trees
    // We focus on *Data mining*.

protected:
    bool _showlog = false;

    // todo: move the filesystem part to the base FormanGradient class
    string _infile;
    string _file_name;
    string _workspace;

    std::map<implicitS, int> _min2lbl;
    // min -> label. valid label starts from 1, <0: invalid labels. 0: not labeled.

    std::vector<int> _pts_lbls;
    size_t _segs_num = 0;

    void _set_workspace(const string &infile) {
        // todo: move to the Formangradient
        _infile = infile;
        std::filesystem::path absolutePath = std::filesystem::absolute(infile);
        _file_name = absolutePath.stem().string();
        std::filesystem::path parentPath = absolutePath.parent_path();
        _workspace = parentPath.string() + "/";

        if (_showlog) {
            cout << "infile: " << _infile << endl;
            cout << "file name: " << _file_name << endl;
            cout << "workspace: " << _workspace << endl;
        }
    }

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
        // Todo: because base class will be initialized at the beginning
        //  need to double check the set_workspace step.
        _showlog = in_debug;
        _set_workspace(infile);
        if (_init() > 0) {
            cerr << "Check infile!\nformat: xyzv\n";
            exit(1);
        }
    }

    // point clustering by Forman over-segmentation
    int cluster(string &out_ptsfile, string &out_minsfile) {
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
        //_pts_lbls = _label_points_by_grouped_mins(gp_mins);
        _label_points_by_grouped_mins_parallel(gp_mins);
        time.stop();
        cout << "pts labeled by mins computed " << time.getElapsedTime() << " s" << endl;

        if (out_ptsfile.empty()) {
            time.start();
            out_ptsfile = _workspace + _file_name + "_lbl.pts";
            write_results_xyzl_txt(out_ptsfile);
            time.stop();
            cout << "pts labeled by mins written " << time.getElapsedTime() << " s" << endl;
        }

        return 0;
        // skip the min network to increase the speed.

        //        std::vector<int> tmp_lbls = _label_points_by_grouped_mins(gp_mins);
        //        string tmp_outfile = _workspace + _file_name + "_lbl_init.vtk";
        //        cout << "write: " << tmp_outfile << endl;
        //        _output_pts_with_label_vtk(tmp_outfile, tmp_lbls, false);


        // create min2mins
        for (const auto &saddle: criticalS[1]) {
            // edge saddles
            vector<implicitS> minima;
            saddle2minima(saddle, minima);
            implicitS m1, m2;
            m1 = minima[0];
            m2 = minima[1];
            if (m1 == m2) {
                continue;
            } else {
                min2mins[m1].insert(m2);
                min2mins[m2].insert(m1);
            }
        } // create min2mins by saddles

        // output file .off: xyzl
        if (out_minsfile.empty()) {
            out_minsfile = _workspace + _file_name + "_mins.off";
        }
        size_t edges_num = 0;
        for (const auto &m: min2mins) {
            edges_num = edges_num + m.second.size();
        }
        ofstream ofs(out_minsfile);
        ofs << fixed << setprecision(3);
        ofs << "# OFF\n";
        ofs << all_mins.size() << " " << edges_num << " 0\n";
        // points
        for (const auto &m: all_mins) {
            vector<int> tmp_vid = m.getConstVertices();
            vector<float> tmp_coord = sc.getVertex(tmp_vid[0]).getCoordinates();
            int tmp_lbl = min2lbl[m];
            ofs << tmp_coord[0] << " " << tmp_coord[1] << " " << tmp_coord[2] << " " << tmp_lbl << endl;
        }
        // edges
        for (const auto &m: min2mins) {
            int mid1 = min2lbl[m.first] - 1; // index = lbl-1
            for (const auto &m2: m.second) {
                int mid2 = min2lbl[m2] - 1;
                ofs << "2 " << mid1 << " " << mid2 << endl;
            }
        }
        ofs.close();


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
    int write_results_xyzl_txt(std::string &outfile = (string &) "") {
        if (outfile.empty()) {
            outfile = _workspace + _file_name + "_lbl.pts";
        }
        cout << "write: " << outfile << endl;
        _output_pts_with_label_pts(outfile, _pts_lbls, false);

        return 0;
    }

    int write_results_xyzl_vtk(string &outfile = (string &) "") {
        if (outfile.empty()) {
            outfile = _workspace + _file_name + "_lbl.vtk";
        }
        cout << "write: " << outfile << endl;
        _output_pts_with_label_vtk(outfile, _pts_lbls, false);

        return 0;
    }

    string output_labeled_pts(const string &outfile);


    // Other functions

    // todo: may independent of tts_aux.h. rename the function
    //      make it as the inline function in toposegment class?
    //      move out of the protect
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
