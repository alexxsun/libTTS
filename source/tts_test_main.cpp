#include <iostream>

#include "projects/xx_tts.h"

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
    return wood_leaf_separation(infile);
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
        if (argc > 4) { // note: -tts is also counted in argc
            th_p2trunk_distance = stod(argv[3]);
        }
        if (argc > 5) {
            th_search_radius = stod(argv[4]);
        }
        test_tree_extraction(infile, trunk_file, th_p2trunk_distance, th_search_radius);
    }

    return 0;
}
