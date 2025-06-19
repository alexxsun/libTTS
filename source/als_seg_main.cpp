#include <iostream>
#include <algorithm>

#include "../source/projects/als_seg.h"

using namespace std;

class InputParser {
    // https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
public:
    InputParser(int &argc, char **argv) {
        for (int i = 1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }

    /// @author iain
    const std::string &getCmdOption(const std::string &option) const {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
            return *itr;
        }
        static const std::string empty_string("");
        return empty_string;
    }

    /// @author iain
    bool cmdOptionExists(const std::string &option) const {
        return std::find(this->tokens.begin(), this->tokens.end(), option)
               != this->tokens.end();
    }

private:
    std::vector<std::string> tokens;
};

int main(int argc, char **argv) {
    // input parameters
    InputParser input(argc, argv);
    if (input.cmdOptionExists("-h")) {
        // Do stuff
        cout << "Usage: " << argv[0] << " xyz.pts -ht -r -h -Fp -pH -pI -ds" << endl;
        cout << "-f xyz.pts or xyz.pts: input vegetation points\n";
        cout << "-ht: type of height function.\n"
                "\t b: height function (clear trunks), i: inverse height function (clear treetops)";
        cout << "Optional parameters\n";
        cout << "-r: radius size. unit: m. default: auto-set.\n";
        cout << "-h: height rate [0,1]. default: 0.5.\n";
        cout << "-Fp: forest-level persistence value [0,1]. default: 0.02.\n";
        cout << "-pH, -pI: tree cluster level persistence value [0,1].\n"
                "\t default values based on height function (ht):\n"
                "\t ht, pH, pI\n"
                "\t b, 0.01, 0.7\n"
                "\t i, 0.7 , 0.01\n";
        cout << "-ds: data shrinking (remove tiny segments based on size). default: > avg_size.";
        cout << "Example\n";
        cout << "Usage: " << argv[0] << " -f xyz.pts -ht b " << "\tclear stems." << endl;
        cout << "Usage: " << argv[0] << " -f xyz.pts -ht t " << "\tclear treetops." << endl;
        cout << "Usage: " << argv[0] << " -f xyz.pts -ht t -ds none" << "\tclear treetops. no postprocess" << endl;
        return 1;
    }

    // required parameters
    string pts_file = argv[1]; // input forest points. xyz.
    const std::string &filename = input.getCmdOption("-f");
    if (!filename.empty()) {
        pts_file = filename;
    } else {
        cout << "require the input file\n";
        return 1;
    }
    cout << "Input pts file: " << pts_file << endl;
    // b: true (bottom up, height function), t: false, otherwise: false
    bool th_height_funtype = true;
    const std::string &height_type = input.getCmdOption("-ht");
    if (!height_type.empty()) {
        th_height_funtype = height_type == "b";
    } else {
        cout << "need to set height function type\n";
    }
    cout << "height function: " << th_height_funtype << endl;

    // optional parameters.
    // radius: <0 auto-set. unit: m. default auto-set
    double th_radius = -1;
    const std::string &tmp = input.getCmdOption("-r");
    if (!tmp.empty()) {
        th_radius = stod(tmp);
    }
    cout << "radius: " << th_radius << endl;

    // height: <0 auto-set [0,1] default 0.5.
    double th_height = 0.5;
    const std::string &tmp_h = input.getCmdOption("-h");
    if (!tmp_h.empty()) {
        th_height = stod(tmp_h);
    }
    cout << "height: " << th_height << endl;

    // persistence values for forest segmentation. <0 auto-set [0,1] default: 0.02.
    double th_forest_pers = 0.02;//stod(argv[3]);
    const std::string &tmp_fp = input.getCmdOption("-Fp");
    if (!tmp_fp.empty()) {
        th_forest_pers = stod(tmp_fp);
    }
    cout << "forest pers: " << th_forest_pers << endl;
    // persistence values for cluster segmentation.
    // <0 auto-set [0,1] default: 0.7 or 0.01 based on type
    double th_pers_H, th_pers_I;
    if (th_height_funtype) {
        th_pers_H = 0.01;
        th_pers_I = 0.7;
    } else {
        th_pers_H = 0.7;
        th_pers_I = 0.01;
    }
    const std::string &tmp_p1 = input.getCmdOption("-pH");
    if (!tmp_p1.empty()) {
        th_pers_H = stod(tmp_p1);
    }
    const std::string &tmp_p2 = input.getCmdOption("-pI");
    if (!tmp_p2.empty()) {
        th_pers_I = stod(tmp_p2);
    }
    cout << "pers H, I: " << th_pers_H << " " << th_pers_I << endl;

    // how to post-process based on single tree pts numbers
    string ds = "newfor";  // how to post-process based on single tree pts numbers
    const std::string &tmp_ds = input.getCmdOption("-ds");
    if (!tmp_ds.empty()) {
        ds = tmp_ds;
    }
    cout << "Data post-processing mode: " << ds << endl;


    //
    cout << "\nProcess data\n";
    process(pts_file, th_radius, th_forest_pers, th_height, th_height_funtype, th_pers_H, th_pers_I, ds);

    cout << "Done\n";
    return 0;
}
