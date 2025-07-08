#include <iostream>

#include "projects/xx_tts.h"
#include "projects/als_seg.h"

using namespace std;


int test_alpha_shape_segment(const string &infile, const double &alpha_sq_value, string &out_tocfile) {
    std::cout << "\nTest alpha shape generation and segmentation\n\n";
    return alpha_shape_segment(infile, alpha_sq_value, out_tocfile);
}

int test_alpha_shape_generation(const string &infile, const double &alpha_sq_value, string &outfile) {
    std::cout << "\nTest alpha shape generation ONLY\n\n";
    return alpha_shape_generation(infile, alpha_sq_value, outfile);
}

int test_wls(const string &infile) {
    std::cout << "\nTest wood-leaf separation. part 1\n\n";
    std::string outfile1, outfile2;
    return get_oversegments(infile, outfile1, outfile2);
}


int test_tree_extraction(const string &vg_meshfile, const string &trunk_file,
                         const double &th_p2trunk_distance,
                         const double &th_search_radius) {
    // with detected tree locations
    std::cout << "\nTest tree extraction.\n\n";
    cout << "th_p2trunk_distance=" << th_p2trunk_distance << "\n";
    cout << "th_search_radius=" << th_search_radius << "\n";
    string outfile;
    return extract_single_trees(vg_meshfile, trunk_file, outfile, th_p2trunk_distance, th_search_radius);
}

int test_als_segment(const string &infile, double th_radius, double th_forest_pers, double th_height,
                     bool th_height_funtype,
                     double th_pers_H, double th_pers_I, const string &ds) {
    return als_segment(infile, th_radius, th_forest_pers, th_height, th_height_funtype,
                       th_pers_H, th_pers_I, ds);
}

int main(int argc, char *argv[]) {
    string last = argv[argc - 1];
    if (last == "-as") {
        string infile = argv[1];
        double alpha_sq = stod(argv[2]);
        string outfile;
        test_alpha_shape_segment(infile, alpha_sq, outfile);
    }

    if (last == "-as2") {
        string infile = argv[1];
        double alpha_sq = stod(argv[2]);
        string outfile; // .off file, or ply file
        test_alpha_shape_generation(infile, alpha_sq, outfile);
    }

    if (last == "-wls") {
        string infile = argv[1];
        test_wls(infile);
    }

    if (last == "-tts") {
        string infile = argv[1]; // .off, .ply
        string trunk_file = argv[2]; // .pts
        double th_p2trunk_distance = 0.2;
        double th_search_radius = 0.25;
        if (argc > 4) {
            // note: -tts is also counted in argc
            th_p2trunk_distance = stod(argv[3]);
        }
        if (argc > 5) {
            th_search_radius = stod(argv[4]);
        }
        test_tree_extraction(infile, trunk_file, th_p2trunk_distance, th_search_radius);
    }

    if (last == "-als") {
        string infile = argv[1]; // .pts file
        double th_radius = -1; // <0 auto-set. unit: m. default auto-set
        double th_forest_pers = 0.02; // <0 auto-set [0,1] default: 0.02.
        double th_height = 0.5; // <0 auto-set [0,1] default 0.5.
        bool th_height_funtype = true; // b: true (bottom up, height function), t: false, otherwise: false
        double th_pers_H = 0.01; // persistence values for cluster segmentation.
        double th_pers_I = 0.7; // persistence values for cluster segmentation.
        string ds = "none"; // how to post-process based on single tree pts numbers

        if (argc > 3) {
            th_radius = stod(argv[2]);
        }
        if (argc > 4) {
            th_forest_pers = stod(argv[3]);
        }
        if (argc > 5) {
            th_height = stod(argv[4]);
        }
        if (argc > 6) {
            string height_type = argv[5];
            th_height_funtype = height_type == "false";
        }
        if (argc > 7) {
            th_pers_H = stod(argv[6]);
        }
        if (argc > 8) {
            th_pers_I = stod(argv[7]);
        }
        if (argc > 9) {
            ds = argv[8];
        }

        test_als_segment(infile, th_radius, th_forest_pers, th_height, th_height_funtype,
                         th_pers_H, th_pers_I, ds);
    }

    return 0;
}
